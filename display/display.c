#include "display.h"

#define FIXED_POINT_SHIFT 8
#define fixed_t int32_t
#define FLOAT_TO_FIXED(x) ((fixed_t)((x) * (1 << FIXED_POINT_SHIFT)))
#define FIXED_TO_FLOAT(x) ((float)(x) / (1 << FIXED_POINT_SHIFT))
#define INT32_TO_FIXED(x) ((fixed_t)((x) << FIXED_POINT_SHIFT))
#define FIXED_TO_INT32(x) ((int32_t)(x) >> FIXED_POINT_SHIFT)
#define FIXED_MUL(x, y) (((x) * (y)) >> FIXED_POINT_SHIFT)
#define FIXED_DIV(x, y) (((x) << FIXED_POINT_SHIFT) / (y))

#define LCD_COLUMNS 400
#define LCD_ROWS 240
#define LCD_STRIDE 52
#define LCD_STRIDE_32 (LCD_STRIDE / 4)

const float xFactor = sqrt(3.f/2.f);
const float yFactor = 1.f/sqrt(2.f);

uint8_t* textureData = NULL;

void loadFile(PlaydateAPI* pd)
{
    pd->system->logToConsole("opening file");
    SDFile* file = pd->file->open("./images/texture.data", kFileRead);
    if(file == NULL)
    {
        pd->system->logToConsole("Error open %s", pd->file->geterr());
        return;
    }
    textureData = pd->system->realloc(textureData, 256*256);
    if(textureData == NULL)
    {
        pd->system->logToConsole("Error alloc");
        return;
    }

    if(pd->file->read(file, textureData, 256*256) < 0)
    {
        pd->system->logToConsole("Error read");
        return;
    }
    pd->system->logToConsole("READ OK");
}


static inline int32_t int32_min(int32_t a, int32_t b)
{ 
	return b + ((a-b) & (a-b)>>31);
}

static inline int32_t int32_max(int32_t a, int32_t b)
{ 
	return a - ((a-b) & (a-b)>>31);
}

static inline uint32_t swap(uint32_t n)
{
#if TARGET_PLAYDATE
	//return __REV(n);
	uint32_t result;
	
	__asm volatile ("rev %0, %1" : "=l" (result) : "l" (n));
	return(result);
#else
	return ((n & 0xff000000) >> 24) | ((n & 0xff0000) >> 8) | ((n & 0xff00) << 8) | (n << 24);
#endif
}

static inline void
_drawMaskPattern(uint32_t* p, uint32_t mask, uint32_t color)
{
	if ( mask == 0xffffffff )
		*p = color;
	else
		*p = (*p & ~mask) | (color & mask);
}

static void
drawFragment(uint32_t* row, int x1, int x2, uint32_t color)
{
	if ( x2 < 0 || x1 >= LCD_COLUMNS )
		return;
	
	if ( x1 < 0 )
		x1 = 0;
	
	if ( x2 > LCD_COLUMNS )
		x2 = LCD_COLUMNS;
	
	if ( x1 > x2 )
		return;
	
	// Operate on 32 bits at a time
	
	int startbit = x1 % 32;
	uint32_t startmask = swap((1 << (32 - startbit)) - 1);
	int endbit = x2 % 32;
	uint32_t endmask = swap(((1 << endbit) - 1) << (32 - endbit));
	
	int col = x1 / 32;
	uint32_t* p = row + col;

	if ( col == x2 / 32 )
	{
		uint32_t mask = 0;
		
		if ( startbit > 0 && endbit > 0 )
			mask = startmask & endmask;
		else if ( startbit > 0 )
			mask = startmask;
		else if ( endbit > 0 )
			mask = endmask;
		
		_drawMaskPattern(p, mask, color);
	}
	else
	{
		int x = x1;
		
		if ( startbit > 0 )
		{
			_drawMaskPattern(p++, startmask, color);
			x += (32 - startbit);
		}
		
		while ( x + 32 <= x2 )
		{
			_drawMaskPattern(p++, 0xffffffff, color);
			x += 32;
		}
		
		if ( endbit > 0 )
			_drawMaskPattern(p, endmask, color);
	}
}

float min2(float a, float b)
{
	return (a < b) ? a : b;
}

float max2(float a, float b)
{
	return (a > b) ? a : b;
}

float min3(float a, float b, float c)
{
	return min2(min2(a, b), c);
}

float max3(float a, float b, float c)
{
	return max2(max2(a, b), c);
}

float orientation(float x1, float y1, float x2, float y2, float x3, float y3)
{
	return (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1);
}

Scene* scene = NULL;

void initDisplay(PlaydateAPI* pd)
{
	loadFile(pd);
  if(scene)
  {
    freeScene(&scene);
    scene = NULL;
  }
  scene = createScene(8, 8);

  scene->p3d[0].x = 1.f;  scene->p3d[0].y = -1.f;  scene->p3d[0].z = -1.f;
	scene->p3d[1].x = 1.f;  scene->p3d[1].y = 1.f;  scene->p3d[1].z = -1.f;
	scene->p3d[2].x = -1.f;  scene->p3d[2].y = 1.f;  scene->p3d[2].z = -1.f;
	scene->p3d[3].x = -1.f;  scene->p3d[3].y = -1.f;  scene->p3d[3].z = -1.f;
	scene->p3d[4].x = 1.f;  scene->p3d[4].y = -1.f;  scene->p3d[4].z = 1.f;
	scene->p3d[5].x = 1.f;  scene->p3d[5].y = 1.f;  scene->p3d[5].z = 1.f;
	scene->p3d[6].x = -1.f;  scene->p3d[6].y = 1.f;  scene->p3d[6].z = 1.f;
	scene->p3d[7].x = -1.f;  scene->p3d[7].y = -1.f;  scene->p3d[7].z = 1.f;

	/*scene->points[0].x = 100;  scene->points[0].y = 80; 
  scene->points[1].x = 100;  scene->points[1].y = 180;
  scene->points[2].x = 200;  scene->points[2].y = 100;
  scene->points[3].x = 200;  scene->points[3].y = 200;
	scene->points[4].x = 300;  scene->points[4].y = 80;
  scene->points[5].x = 300;  scene->points[5].y = 180;*/

	scene->uvs[0].u = 0.0f;  scene->uvs[0].v = 0.5f;
  scene->uvs[1].u = 0.5f;  scene->uvs[1].v = 0.5f;
  scene->uvs[2].u = 1.0f;  scene->uvs[2].v = 0.5f;
  scene->uvs[3].u = 0.5f;  scene->uvs[3].v = 0.5f;
	scene->uvs[4].u = 0.0f;  scene->uvs[4].v = 1.0f;
  scene->uvs[5].u = 0.5f;  scene->uvs[5].v = 1.0f;
  scene->uvs[6].u = 1.0f;  scene->uvs[6].v = 1.0f;
  scene->uvs[7].u = 0.5f;  scene->uvs[7].v = 1.0f;

  /*scene->faces[0].points[0] = 0; scene->faces[0].points[1] = 1; scene->faces[0].points[2] = 4;
  scene->faces[1].points[0] = 1; scene->faces[1].points[1] = 5; scene->faces[1].points[2] = 4;
	scene->faces[2].points[0] = 1; scene->faces[2].points[1] = 2; scene->faces[2].points[2] = 5;
  scene->faces[3].points[0] = 2; scene->faces[3].points[1] = 6; scene->faces[3].points[2] = 5;
	scene->faces[4].points[0] = 2; scene->faces[4].points[1] = 3; scene->faces[4].points[2] = 6;
  scene->faces[5].points[0] = 3; scene->faces[5].points[1] = 7; scene->faces[5].points[2] = 6;
	scene->faces[6].points[0] = 3; scene->faces[6].points[1] = 0; scene->faces[6].points[2] = 7;
  scene->faces[7].points[0] = 0; scene->faces[7].points[1] = 4; scene->faces[7].points[2] = 7;*/
	scene->faces[0].points[0] = 0; scene->faces[0].points[2] = 1; scene->faces[0].points[1] = 4;
  scene->faces[1].points[0] = 1; scene->faces[1].points[2] = 5; scene->faces[1].points[1] = 4;
	scene->faces[2].points[0] = 1; scene->faces[2].points[2] = 2; scene->faces[2].points[1] = 5;
  scene->faces[3].points[0] = 2; scene->faces[3].points[2] = 6; scene->faces[3].points[1] = 5;
	scene->faces[4].points[0] = 2; scene->faces[4].points[2] = 3; scene->faces[4].points[1] = 6;
  scene->faces[5].points[0] = 3; scene->faces[5].points[2] = 7; scene->faces[5].points[1] = 6;
	scene->faces[6].points[0] = 3; scene->faces[6].points[2] = 0; scene->faces[6].points[1] = 7;
  scene->faces[7].points[0] = 0; scene->faces[7].points[2] = 4; scene->faces[7].points[1] = 7;
}

Scene* createScene(size_t pointCount, size_t faceCount)
{
  Scene* scene = (Scene*)malloc(sizeof(Scene));
	scene->p3d = (Point3D*)malloc(sizeof(Point3D) * pointCount);
  scene->points = (Point2D*)malloc(sizeof(Point2D) * pointCount);
  scene->uvs = (UV*)malloc(sizeof(UV) * pointCount);
  scene->faces = (Face*)malloc(sizeof(Face) * faceCount);
  scene->pointCount = pointCount;
  scene->faceCount = faceCount;

  return scene;
}

void freeScene(Scene** scene)
{
  if(*scene)
  {
    free((*scene)->points);
    free((*scene)->uvs);
    free((*scene)->faces);
    free(*scene);
    *scene = NULL;
  }
}

void drawFace(uint8_t* frameBuffer, Face* face, PlaydateAPI* pd)
{
  uint32_t* row = (uint32_t*)frameBuffer;
	float x1 = scene->points[face->points[0]].x;
	float y1 = scene->points[face->points[0]].y;
	float x2 = scene->points[face->points[1]].x;
	float y2 = scene->points[face->points[1]].y;
	float x3 = scene->points[face->points[2]].x;
	float y3 = scene->points[face->points[2]].y;

	float u1 = scene->uvs[face->points[0]].u;
	float v1 = scene->uvs[face->points[0]].v;
	float u2 = scene->uvs[face->points[1]].u;
	float v2 = scene->uvs[face->points[1]].v;
	float u3 = scene->uvs[face->points[2]].u;
	float v3 = scene->uvs[face->points[2]].v;
  /*for(int y = y1; y <= y2; y++)
  {
    drawFragment(row + y * LCD_STRIDE_32, (int)x1, (int)x2, 0xffffffff);
  }*/

	// Compute bounding boxes

	float minX = min3(x1, x2, x3);
	float maxX = max3(x1, x2, x3);
	float minY = min3(y1, y2, y3);
	float maxY = max3(y1, y2, y3);

	// Clip against screen bounds;

	minX = max2(minX, 0);
	minY = max2(minY, 0);
	maxX = min2(maxX, LCD_COLUMNS - 1);
	maxY = min2(maxY, LCD_COLUMNS - 1);

  // Rasterize

  // Triangle setup
  float a12 = y1 - y2, b12 = x2 - x1;
  float a23 = y2 - y3, b23 = x3 - x2;
  float a31 = y3 - y1, b31 = x1 - x3;

  float w1_row = orientation(x2, y2, x3, y3, minX, minY);
  float w2_row = orientation(x3, y3, x1, y1, minX, minY);
  float w3_row = orientation(x1, y1, x2, y2, minX, minY);

	float w1 = 0.f;
	float w2 = 0.f;
	float w3 = 0.f;

	for(size_t y = minY + 0.5f; y <= maxY + 0.5f; y++, w1_row += b23, w2_row += b31, w3_row += b12)
	{
		w1 = w1_row;
		w2 = w2_row;
		w3 = w3_row;

		for(size_t x = minX + 0.5f; x <= maxX + 0.5f; x++, w1 += a23, w2 += a31, w3 += a12)
		{
			if(w1 >= 0.f && w2 >= 0.f && w3 >= 0.f)
			{
				float u = (w1 * u1 + w2 * u2 + w3 * u3)/(w1+w2+w3);
				float v = (w1 * v1 + w2 * v2 + w3 * v3)/(w1+w2+w3);
				uint8_t _col = u * 255;
				uint8_t _row = v * 255;
				uint8_t val = textureData[_col + _row * 256];
				drawFragment(row + y * LCD_STRIDE_32, (int)x, (int)x + 1, val == 0 ? 0x00000000 : 0xffffffff);
			}
		}
	}

}


void draw(uint8_t* frameBuffer, PlaydateAPI* pd)
{
  //float angle = 45.f/180.f * M_PI;
	float angle = pd->system->getCrankAngle()/180.f*M_PI;

  if(scene)
  {
		for(size_t i = 0; i < scene->pointCount; i++)
		{
			const float factor = 50.f;
			const float dX = 200.f;
			const float dY = 120.f;
			scene->points[i].x = xFactor * (scene->p3d[i].x * cos(angle) - scene->p3d[i].y * sin(angle));
			scene->points[i].x = factor * scene->points[i].x + dX;
			scene->points[i].y = yFactor * (scene->p3d[i].x * sin(angle) + scene->p3d[i].y * cos(angle)) + scene->p3d[i].z;
			scene->points[i].y = factor * scene->points[i].y + dY;
		}
		pd->system->logToConsole("Angle %f -> (%f, %f)", angle, scene->points[0].x, scene->points[0].y);
    for(size_t i = 0; i < scene->faceCount; i++)
    {
      drawFace(frameBuffer, &scene->faces[i], pd);
    }
  }
	//drawFragment((uint32_t*)frameBuffer + 120 * LCD_STRIDE_32, 0, LCD_COLUMNS, 0xff00ff00);
}
