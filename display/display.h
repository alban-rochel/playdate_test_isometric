#pragma once

#include "pd_api.h"

typedef struct
{
  float x;
  float y;
  float z;
} Point3D;

typedef struct
{
  float x;
  float y;
} Point2D;

typedef struct
{
  float u;
  float v;
} UV;

typedef struct{

  size_t points[3];
  size_t uvs[3];
} Face;

typedef struct{

  size_t points[4];
  size_t uvs[4];
} Quad;

typedef struct
{
  Point3D* p3d;
  Point2D* points;
  UV* uvs;
  Face* faces;
  Quad* quads;
  size_t pointCount;
  size_t faceCount;
  size_t quadCount;
} Scene;

void initDisplay(PlaydateAPI* pd);

Scene* createScene(size_t pointCount, size_t faceCount, size_t quadCount, size_t uvCount);
void freeScene(Scene** scene);

void drawFace(uint8_t* frameBuffer, Face* face, PlaydateAPI* pd);

void draw(uint8_t* frameBuffer, PlaydateAPI* pd);
