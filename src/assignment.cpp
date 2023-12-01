#include <iostream>

#include "image.h"
#include "object.h"
#include "scene.h"

using namespace Eigen;
using namespace std;

const int MAX_ITERS = 10000;
const int XRES = 1000;
const int YRES = 1000;

const double epsilon = 1e-4;

/**
 * Helpers
 */

inline int sign(float x) { return (x > 0) ? 1 : -1; }

inline float rad2deg(float angle) { return angle * 180 / M_PI; }

inline float deg2rad(float angle) { return angle * M_PI / 180; }

inline double superquadricIO(Vector3d point, double e, double n) {
  double x2 = pow(point[0], 2);
  double y2 = pow(point[1], 2);
  double z2 = pow(point[2], 2);
  double S =
      -1 + pow(z2, 1.0 / n) + pow(pow(x2, 1.0 / e) + pow(y2, 1.0 / e), e / n);
  return S;
}

inline Vector3d superquadricIOGrad(Vector3d point, double e, double n) {
  double x2 = pow(point[0], 2);
  double y2 = pow(point[1], 2);
  double z2 = pow(point[2], 2);

  double common_term = pow(pow(x2, 1.0 / e) + pow(y2, 1.0 / e), (e / n) - 1.0);

  Vector3d grad = {2 * point[0] * pow(x2, (1.0 / e) - 1.0) * common_term,
                   2 * point[1] * pow(y2, (1.0 / e) - 1.0) * common_term,
                   2 * point[2] * pow(z2, (1.0 / n) - 1.0)};
  return grad / n;
}
/**
 * IOTest Code
 */

bool Superquadric::IOTest(const Vector3d &point) {
  /**
   * PART 1
   * TODO: Implement the IO Test function for a superquadric. Remember to
   *       apply any transformations to the superquadric before testing
   *       against the IO function.
   */
  /* Copy and convert to homogeneous coordinates */
  Vector4d p = {point[0], point[1], point[2], 1};

  /* Transform to body coordinates */
  for (auto it = transforms.rbegin(); it != transforms.rend(); it++) {
    p = (*it)->GetMatrix().inverse() * p;
  }

  /* Convert back from homogeneous coordinates */
  const Vector3d new_point(p[0] / p[3], p[1] / p[3], p[2] / p[3]);

  /* Calculate IO function*/
  double S = superquadricIO(new_point, exp0, exp1);

  /* Return true if IO function <= 0, false if > 0 */
  return S <= 0;
}

bool Assembly::IOTest(const Vector3d &point) {
  /**
   * PART 1
   * TODO: Implement the IO Test function for an assembly (recursively call
   *       IOTest on the children). Make sure to apply any transformations
   *       to the assembly before calling IOTest on the children.
   */
  /* Copy and convert to homogeneous coordinates */
  Vector4d p = {point[0], point[1], point[2], 1};

  /* Transform to body coordinates */
  for (auto it = transforms.rbegin(); it != transforms.rend(); it++) {
    p = (*it)->GetMatrix().inverse() * p;
  }

  /* Convert back from homogeneous coordinates */
  const Vector3d new_point(p[0] / p[3], p[1] / p[3], p[2] / p[3]);

  /* Call IOTest on children */
  for (auto &child : children) {
    if (child->IOTest(new_point))
      return true;
  }
  return false;
}

/**
 * Closest Intersection Code
 */

pair<double, Intersection> Superquadric::ClosestIntersection(const Ray &ray) {
  /**
   * PART 1
   * TODO: Implement a ray-superquadric intersection using Newton's method.
   *       Make sure to apply any transformations to the superquadric before
   *       performing Newton's method.
   */
  /* Transform ray to body space */
  Ray incoming = ray;
  for (auto it = transforms.rbegin(); it != transforms.rend(); it++) {
    incoming.Transform((*it)->GetMatrix().inverse());
  }

  /* Initialize return pair */
  pair<double, Intersection> closest = make_pair(INFINITY, Intersection());

  /* Initial guess using bounding sphere */
  Vector3d av = incoming.direction;
  Vector3d bv = incoming.origin;

  double a = av.dot(av);
  double b = 2 * av.dot(bv);
  double c = bv.dot(bv) - 3;

  double determinant = pow(b, 2) - 4 * a * c;
  if (determinant < 0)
    return closest;

  double common_term = -b - sign(b) * sqrt(determinant);
  double t1 = common_term / (2 * a);
  double t2 = 2 * c / common_term;

  // If we might havea hit...
  if (t1 > 0 && t2 > 0) {
    // Choose closer intersection
    double t = min(t1, t2);

    // Execute Newton's method
    double g = INFINITY;
    double gp = -INFINITY;
    while (abs(g) > epsilon && gp <= 0) {
      g = superquadricIO(incoming.At(t), exp0, exp1);
      gp = incoming.direction.dot(
          superquadricIOGrad(incoming.At(t), exp0, exp1));
      t = t - (g / gp);
    }

    // Check if we're hitting
    g = superquadricIO(incoming.At(t), exp0, exp1);
    if (abs(g) <= epsilon) {
      // Define ray for intersection location, with direction being that
      // of the reflected ray
      Ray intersection;
      intersection.origin = incoming.At(t);
      Vector3d d = incoming.direction;
      Vector3d n = superquadricIOGrad(incoming.At(t), exp0, exp1).normalized();
      intersection.direction = d - 2 * d.dot(n) * n;
      intersection.Normalize();

      // Need to transform back to world coordinates
      for (auto it = transforms.begin(); it != transforms.end(); it++) {
        intersection.Transform((*it)->GetMatrix());
      }

      // Save
      closest.first = t;
      closest.second.location = intersection;
      closest.second.obj = this;
    } // Else we've missed
  }   // Else we've missed

  return closest;
}

pair<double, Intersection> Assembly::ClosestIntersection(const Ray &_ray) {
  /**
   * PART 1
   * TODO: Implement a ray-assembly intersection by recursively finding
   *       intersection with the assembly's children. Make sure to apply any
   *       transformations to the assembly before calling ClosestIntersection
   *       on the children.
   */
  /* Transform ray to body space */
  Ray ray = _ray;
  for (auto it = transforms.rbegin(); it != transforms.rend(); it++) {
    ray.Transform((*it)->GetMatrix().inverse());
  }

  /* Find children with closest intersection */
  double mint = INFINITY;
  Intersection minIntersection = Intersection();
  for (auto &child : children) {
    auto closest = child->ClosestIntersection(ray);
    if (closest.first < mint) {
      mint = closest.first;
      minIntersection = closest.second;
    }
  }

  // Transform intersection to current-body coordinates
  for (auto it = transforms.begin(); it != transforms.end(); it++) {
    minIntersection.location.Transform((*it)->GetMatrix());
  }

  return make_pair(mint, minIntersection);
}

/**
 * Raytracing Code
 */

void Scene::Raytrace() {
  // Initialize image for drawing
  Image img = Image(XRES, YRES);

  // Get scene camera
  Camera camera = GetCamera();
  Vector3d camera_pos = -camera.translate.GetDelta();

  // Find width and height of screen (back of frustum) in world space
  double height = 2 * camera.frustum.near * tan(deg2rad(camera.frustum.fov / 2));
  double width = camera.frustum.aspect_ratio * height;

  // Find basis bectors while transforming them to world space
  Vector4d _e[] = {{0, 0, -1, 1}, {1, 0, 0, 1}, {0, 1, 0, 1}};
  for (int i = 0; i < 3; i++) {
    _e[i] = camera.rotate.GetMatrix().inverse() * _e[i];
  }
  Vector3d e1 = _e[0].head(3) / _e[0][3];
  Vector3d e2 = _e[1].head(3) / _e[1][3];
  Vector3d e3 = _e[2].head(3) / _e[2][3];

  for (int i = 0; i < XRES; i++) {
    printf("i = %d\n", i);
    for (int j = 0; j < YRES; j++) {
      /**
       * PART 2
       * TODO: Implement raytracing using the code from the first part
       *       of the assignment. Set the correct color for each pixel
       *       here.
       */
      // Initialize color of pixel
      Vector3f color = Vector3f::Zero();

      // Find ray
      Ray ray;
      ray.origin = camera_pos;

      double xi = (i - XRES / 2) * (width / XRES);
      double yj = (j - YRES / 2) * (height / YRES);
      ray.direction = camera.frustum.near * e1 + xi * e2 + yj * e3;

      // "Send it out"
      auto closest = ClosestIntersection(ray);
      if (closest.first != INFINITY) {
        // Define variables about interseciton...
        Superquadric *obj = closest.second.obj; // intersecting object
        Material material = obj->GetMaterial(); // its material
        Vector3d intersection_pos = closest.second.location.origin;

        // Initialize color accumulators
        Vector3f diffuse_sum = Vector3f::Zero();
        Vector3f specular_sum = Vector3f::Zero();

        // We have a ray-superquadric intersection
        // For each light, see if its incoming ray isn't obstructed
        for (Light light : lights) {
          // Get light position as Vector3d
          Vector3d light_pos = light.position.head(3) / light.position[3];

          // "Send out" a light ray
          Ray light_ray;
          light_ray.origin = light_pos;
          light_ray.direction = intersection_pos - light_pos;
          auto closest2 = ClosestIntersection(light_ray);

          // If the closest intersection from the light is the same
          // intersection...
          if ((closest.second.location.origin - closest2.second.location.origin)
                  .norm() < 2*epsilon) {
            // Calculate lighting model
            Vector3d light_vec = light_pos - intersection_pos;
            Vector3d camera_vec = camera_pos - intersection_pos;

            // Find normal
            Vector3d normal = obj->GetNormal(intersection_pos);

            // Find attenuation coefficient
            double distance = light_vec.norm();
            double coef = 1.0 / (1 + light.attenuation * pow(distance, 2));

            // Diffuse component
            diffuse_sum += coef * light.color.ToVector() * max(0.0, normal.dot(light_vec));

            // Specular component
            specular_sum += coef * light.color.ToVector() *
                            pow(max(0.0, normal.dot((camera_vec.normalized() +
                                                     light_vec.normalized())
                                                        .normalized())),
                                material.shininess);
          }
        }
        color = material.ambient.ToVector() +
                material.diffuse.ToVector().cwiseProduct(diffuse_sum) +
                material.specular.ToVector().cwiseProduct(specular_sum);
      }
      img.SetPixel(i, j, color);
    }
  }

  // Outputs the image.
  if (!img.SaveImage("rt.png")) {
    cerr << "Error: couldn't save PNG image" << std::endl;
  } else {
    cout << "Done!\n";
  }
}
