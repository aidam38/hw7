#include <iostream>

#include "image.h"
#include "object.h"
#include "scene.h"

using namespace Eigen;
using namespace std;

const int MAX_ITERS = 10000;
const int XRES = 500;
const int YRES = 500;

/**
 * Helpers
 */

inline int sign(float x) { return (x > 0) ? 1 : -1; }

inline double superquadricIO(Vector3d point, double e, double n) {
  double x2 = pow(point[0], 2);
  double y2 = pow(point[1], 2);
  double z2 = pow(point[2], 2);
  double S =
      -1 + pow(z2, 1.0 / n) + pow(pow(x2, 1.0 / e) + pow(y2, 1.0 / e), e / n);
  return S;
}

inline superquadricIOGrad(Vector3d point, double e, double n) {
  double x2 = pow(point[0], 2);
  double y2 = pow(point[1], 2);
  double z2 = pow(point[2], 2);

  double common_term = pow(pow(x2, 1.0 / e) + pow(y2, 1.0 / e), e / n - 1);

  Vector3d grad = {2 * point[0] * pow(x2, 1.0 / e - 1) * common_term,
                   2 * point[1] * pow(y2, 1.0 / e - 1) * common_term,
                   2 * point[2] * pow(z2, 1.0 / e - 1)};
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
  for (auto &t : transforms) {
    p = t->GetMatrix().inverse() * p;
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
  for (auto &t : transforms) {
    p = t->GetMatrix().inverse() * p;
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
  for (auto &t : transforms) {
    ray.Transform(t->GetMatrix().inverse());
  }

  /* Initial guess using bounding sphere */
  Vector3d av = ray.direction;
  Vector3d bv = ray.origin;

  double a = av.dot(av);
  double b = 2 * av.dot(bv);
  double c = bv.dot(bv) - 3;

  double common_term = -b - sign(b) * sqrt(pow(b, 2) - 4 * a * c);
  double t1 = common_term / (2 * a);
  double t2 = 2 * c / common_term;

  /* Recover result */
  pair<double, Intersection> closest = make_pair(INFINITY, Intersection());

  if (t1 > 0 && t2 > 0) {
    // Choose closer intersection
    double t = min(t1, t2);

    // Execute Newton's method
    while (true) {
      t = t - superquadricIO(ray.At(t), exp0, exp1) /
                  superquadricIOGrad(ray.At(t), exp0, exp1);
    }

  } // Else we've missed completely

  return closest;
}

pair<double, Intersection> Assembly::ClosestIntersection(const Ray &ray) {
  /**
   * PART 1
   * TODO: Implement a ray-assembly intersection by recursively finding
   *       intersection with the assembly's children. Make sure to apply any
   *       transformations to the assembly before calling ClosestIntersection
   *       on the children.
   */
  pair<double, Intersection> closest = make_pair(INFINITY, Intersection());
  return closest;
}

/**
 * Raytracing Code
 */

void Scene::Raytrace() {
  Image img = Image(XRES, YRES);

  for (int i = 0; i < XRES; i++) {
    for (int j = 0; j < YRES; j++) {
      /**
       * PART 2
       * TODO: Implement raytracing using the code from the first part
       *       of the assignment. Set the correct color for each pixel
       *       here.
       */
      img.SetPixel(i, j, Vector3f::Ones());
    }
  }

  // Outputs the image.
  if (!img.SaveImage("rt.png")) {
    cerr << "Error: couldn't save PNG image" << std::endl;
  } else {
    cout << "Done!\n";
  }
}
