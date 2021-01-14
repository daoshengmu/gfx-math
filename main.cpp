//
//  main.cpp
//  gfx-math
//
//  Created by Daosheng Mu on 9/26/20.
//  Copyright © 2020 Daosheng Mu. All rights reserved.
//

#include <iostream>
#include "Vector3d.h"
#include "Matrix4x4.h"
#include "Quaternion.h"

using namespace gfx_math;

int main(int argc, const char * argv[]) {
  std::cout << "------Test Start-----" << std::endl;
  std::cout << std::endl;
  
  Vector3D<float> a(1, 0, 0);
  
  float res = a.DotProduct(Vector3D<float>(1, 0, 0));
  assert(res == 1);
  
  Matrix4x4<float> matA;
  Matrix4x4<float> matB;
  matA = matB;
  matA.Identity();
  Matrix4x4<float> matC = matA * matB;
  matC = Matrix4x4<float>::LookAtMatrix(Vector3D<float>(0, 0.5, 0.5),
                                        Vector3D<float>(0,0.2,-1),
                                        Vector3D<float>(0,1,0));
  
  Matrix4x4<float> matD = Matrix4x4<float>::TargetTo(Vector3D<float>(0, 0.5, 0.5),                                                                                           Vector3D<float>(0,0.2,-1),
                                                     Vector3D<float>(0,1,0));
  assert(matC != matD); // TODO: We haven't known how to test it.
  
  Matrix4x4<float> matE = Matrix4x4<float>::FPSView(Vector3D<float>(0, 0.5, 0.5), 0.1, 0.2);
  
  Matrix4x4<float> matF = Matrix4x4<float>::Arcball(Vector3D<float>(0, 0.5, 0.5), Quaternion<float>(1,1,0,1), Vector3D<float>(1.0, 0.5, 0.5));
  assert(matE != matF); // TODO: We haven't known how to test it.
  
  Vector3Df ray(1, 0, 0);
  Vector3Df p0(-2, 0, 0);
  Vector3Df planeNrm(-1, 0, 0);
  float planeDist = 10;
  Vector3Df interPoint;
  bool result = RayIntersectWithPlane(ray, p0, planeNrm, planeDist, interPoint);
  assert(result);
  
  res = PointIntersectWithPlane(interPoint, planeNrm, planeDist);
  assert(res == 0);
  
  float dist = PointDistanceWithPlane(interPoint + Vector3Df(1, 0, 0), planeNrm, planeDist);
  assert(dist == 1);
  
  std::vector<Vector3Df> interPointList;
  result = RayIntersectWithSphere(ray, Vector3Df(-10, 0, 0), Vector3Df(0,0,0), 5.0f, interPointList);
  assert(result);
  assert(interPointList.size() == 2);
  
  interPointList.clear();
  result = RayIntersectWithSphere(ray, p0, Vector3Df(0,0,0), 5.0f, interPointList);
  assert(result);
  assert(interPointList.size() == 1);
  
  Vector3Df p;
  result = RayIntersectWithTriangle(Vector3Df(0,0,-1), p0, Vector3Df(-10,-1,-5), Vector3Df(10,-1,-5), Vector3Df(5,1,-5), p);
  assert(result);
  assert(p == Vector3Df(-2, 0, -5));
  
  result = RayIntersectWithAABBox(ray, p0, Vector3Df(2, -2, -2), Vector3Df(6, 2, 2), p);
  assert(result);
  assert(p == Vector3Df(2, 0, 0));
  
  result = SphereIntersectWithAABBox(Vector3Df(0,0,0), 5.0f, Vector3Df(2, -2, -2), Vector3Df(6, 2, 2));
  assert(result);
  
  result = AABBoxIntersectWithAABBox(Vector3Df(0, 0, -2), Vector3Df(4, 4, 2), Vector3Df(2, -2, -2), Vector3Df(6, 2, 2));
  assert(result);

  p = ReflectionVector(Vector3Df(0.707, 0.707, 0), Vector3Df(0, 1, 0));
  assert(p == Vector3Df(-0.707, 0.707, 0));
  
  std::cout << "------Test End-------" << std::endl;
  std::cout << std::endl;
  
  return 0;
}
