#pragma once

#include "pd_api.h"

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

typedef struct
{
  Point2D* points;
  UV* uvs;
  Face* faces;
  size_t pointCount;
  size_t faceCount;
} Scene;

void initDisplay(PlaydateAPI* pd);

Scene* createScene(size_t pointCount, size_t faceCount);
void freeScene(Scene** scene);

void drawFace(uint8_t* frameBuffer, Face* face, PlaydateAPI* pd);

void draw(uint8_t* frameBuffer, PlaydateAPI* pd);
