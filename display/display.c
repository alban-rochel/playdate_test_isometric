#include "display.h"
#include <string.h>

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
uint8_t* __frameBuffer__ = NULL;

void loadFile(PlaydateAPI* pd)
{
    pd->system->logToConsole("opening file");
    SDFile* file = pd->file->open("./images/texture.data", kFileRead);
    if(file == NULL)
    {
        pd->system->logToConsole("Error open %s", pd->file->geterr());
        return;
    }
    textureData = pd->system->realloc(textureData, 128*128);
    if(textureData == NULL)
    {
        pd->system->logToConsole("Error alloc");
        return;
    }

    if(pd->file->read(file, textureData, 128*128) < 0)
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

float min4(float a, float b, float c, float d)
{
	return min2(min2(a, b), min2(c, d));
}

float max4(float a, float b, float c, float d)
{
	return max2(max2(a, b), max2(c, d));
}

float orientation(float x1, float y1, float x2, float y2, float x3, float y3)
{
	return (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1);
}

Scene* scene = NULL;

void initDisplay(PlaydateAPI* pd)
{
	__frameBuffer__ = pd->system->realloc(__frameBuffer__, 400*240);
	loadFile(pd);
  if(scene)
  {
    freeScene(&scene);
    scene = NULL;
  }
  scene = createScene(10, 2, 6/*14*/);

  scene->p3d[0].x = 1.5f;  scene->p3d[0].y = -1.f;  scene->p3d[0].z = -1.f;
	scene->p3d[1].x = 1.5f;  scene->p3d[1].y = 1.f;  scene->p3d[1].z = -1.f;
	scene->p3d[2].x = -1.5f;  scene->p3d[2].y = 1.f;  scene->p3d[2].z = -1.f;
	scene->p3d[3].x = -1.5f;  scene->p3d[3].y = -1.f;  scene->p3d[3].z = -1.f;
	scene->p3d[4].x = 1.5f;  scene->p3d[4].y = -1.f;  scene->p3d[4].z = 1.f;
	scene->p3d[5].x = 1.5f;  scene->p3d[5].y = 1.f;  scene->p3d[5].z = 1.f;
	scene->p3d[6].x = -1.5f;  scene->p3d[6].y = 1.f;  scene->p3d[6].z = 1.f;
	scene->p3d[7].x = -1.5f;  scene->p3d[7].y = -1.f;  scene->p3d[7].z = 1.f;

	scene->p3d[8].x = 0.f;  scene->p3d[8].y = 1.f;  scene->p3d[8].z = -2.f;
	scene->p3d[9].x = 0.f;  scene->p3d[9].y = -1.f;  scene->p3d[9].z = -2.f;

	scene->uvs[0].u = 0.0f * 127.f;  scene->uvs[0].v = 0.5f * 127.f;
  scene->uvs[1].u = 0.5f * 127.f;  scene->uvs[1].v = 0.5f * 127.f;
  scene->uvs[2].u = 1.0f * 127.f;  scene->uvs[2].v = 0.5f * 127.f;
  scene->uvs[3].u = 0.5f * 127.f;  scene->uvs[3].v = 0.5f * 127.f;
	scene->uvs[4].u = 0.0f * 127.f;  scene->uvs[4].v = 1.0f * 127.f;
  scene->uvs[5].u = 0.5f * 127.f;  scene->uvs[5].v = 1.0f * 127.f;
  scene->uvs[6].u = 1.0f * 127.f;  scene->uvs[6].v = 1.0f * 127.f;
  scene->uvs[7].u = 0.5f * 127.f;  scene->uvs[7].v = 1.0f * 127.f;

	scene->uvs[8].u = 0.0f * 127.f;  scene->uvs[8].v = 0.0f * 127.f;
  scene->uvs[9].u = 0.0f * 127.f;  scene->uvs[9].v = 0.0f * 127.f;

	/*scene->faces[0].points[0] = 0; scene->faces[0].points[1] = 4; scene->faces[0].points[2] = 1;
  scene->faces[1].points[0] = 1; scene->faces[1].points[1] = 4; scene->faces[1].points[2] = 5; */
	scene->quads[0].points[0] = 0; scene->quads[0].points[1] = 4; scene->quads[0].points[2] = 5; scene->quads[0].points[3] = 1;

	/*scene->faces[2].points[0] = 1; scene->faces[2].points[1] = 5; scene->faces[2].points[2] = 2;
  scene->faces[3].points[0] = 2; scene->faces[3].points[1] = 5; scene->faces[3].points[2] = 6;*/
	scene->quads[1].points[0] = 1; scene->quads[1].points[1] = 5; scene->quads[1].points[2] = 6; scene->quads[1].points[3] = 2;

	/*scene->faces[4].points[0] = 2; scene->faces[4].points[1] = 6; scene->faces[4].points[2] = 3;
  scene->faces[5].points[0] = 3; scene->faces[5].points[1] = 6; scene->faces[5].points[2] = 7;*/
	scene->quads[2].points[0] = 2; scene->quads[2].points[1] = 6; scene->quads[2].points[2] = 7; scene->quads[2].points[3] = 3;
	
	/*scene->faces[6].points[0] = 3; scene->faces[6].points[1] = 7; scene->faces[6].points[2] = 0;
  scene->faces[7].points[0] = 0; scene->faces[7].points[1] = 7; scene->faces[7].points[2] = 4;*/
	scene->quads[3].points[0] = 3; scene->quads[3].points[1] = 7; scene->quads[3].points[2] = 4; scene->quads[3].points[3] = 0;

	scene->faces[0].points[0] = 1; scene->faces[0].points[1] = 2; scene->faces[0].points[2] = 8;
  scene->faces[1].points[0] = 3; scene->faces[1].points[1] = 0; scene->faces[1].points[2] = 9;

  /*scene->faces[10].points[0] = 0; scene->faces[10].points[1] = 1; scene->faces[10].points[2] = 8;
  scene->faces[11].points[0] = 8; scene->faces[11].points[1] = 9; scene->faces[11].points[2] = 0;*/
	scene->quads[4].points[0] = 1; scene->quads[4].points[1] = 8; scene->quads[4].points[2] = 9; scene->quads[4].points[3] = 0;

  /*scene->faces[12].points[0] = 2; scene->faces[12].points[1] = 3; scene->faces[12].points[2] = 8;
  scene->faces[13].points[0] = 8; scene->faces[13].points[2] = 9; scene->faces[13].points[1] = 3;*/
	scene->quads[5].points[0] = 2; scene->quads[5].points[1] = 3; scene->quads[5].points[2] = 9; scene->quads[5].points[3] = 8;
}

Scene* createScene(size_t pointCount, size_t faceCount, size_t quadCount)
{
  Scene* scene = (Scene*)malloc(sizeof(Scene));
	scene->p3d = (Point3D*)malloc(sizeof(Point3D) * pointCount);
  scene->points = (Point2D*)malloc(sizeof(Point2D) * pointCount);
  scene->uvs = (UV*)malloc(sizeof(UV) * pointCount);
  scene->faces = (Face*)malloc(sizeof(Face) * faceCount);
  scene->quads = (Quad*)malloc(sizeof(Quad) * quadCount);
  scene->pointCount = pointCount;
  scene->faceCount = faceCount;
  scene->quadCount = quadCount;

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
	// Minimize float ops, hoist invariants, use integer math where possible
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

	// Compute bounding box (use ints for loop bounds)
	float minXf = min3(x1, x2, x3);
	float maxXf = max3(x1, x2, x3);
	float minYf = min3(y1, y2, y3);
	float maxYf = max3(y1, y2, y3);

	int minX = (int)max2(minXf, 0);
	int minY = (int)max2(minYf, 0);
	int maxX = (int)min2(maxXf, LCD_COLUMNS - 1);
	int maxY = (int)min2(maxYf, LCD_ROWS - 1);

	// Precompute triangle setup
	float a12 = y1 - y2, b12 = x2 - x1;
	float a23 = y2 - y3, b23 = x3 - x2;
	float a31 = y3 - y1, b31 = x1 - x3;

	float w12_row = orientation(x1, y1, x2, y2, minX + 0.5f, minY + 0.5f);
	float w23_row = orientation(x2, y2, x3, y3, minX + 0.5f, minY + 0.5f);
	float w31_row = orientation(x3, y3, x1, y1, minX + 0.5f, minY + 0.5f);

	float invw123 = 1.f / (w12_row + w23_row + w31_row);

	// Precompute UV deltas for scanline stepping
	float deltaU132 = (a23 * u1 + a31 * u2 + a12 * u3) * invw123;
	float deltaV132 = (a23 * v1 + a31 * v2 + a12 * v3) * invw123;

	for(int y = minY; y <= maxY; ++y,
			w23_row += b23, w31_row += b31, w12_row += b12)
	{
		float w23 = w23_row;
		float w31 = w31_row;
		float w12 = w12_row;

		uint8_t* target = &frameBuffer[y * 400 + minX];

		float u132 = (w23 * u1 + w31 * u2 + w12 * u3) * invw123;
		float v132 = (w23 * v1 + w31 * v2 + w12 * v3) * invw123;

		int done = 0;
		for(int x = minX; x <= maxX; ++x,
				w23 += a23, w31 += a31, w12 += a12,
				++target,
				u132 += deltaU132, v132 += deltaV132)
		{
			if(w31 >= 0.f && w23 >= 0.f && w12 >= 0.f)
			{
				int _col = (int)u132;
				int _row = (int)v132;
				// Clamp UVs to texture size
				//if(_col >= 0 && _col < 128 && _row >= 0 && _row < 128)
					*target = textureData[_col + _row * 128];
				done = 1;
			}
			else if(done)
			{
				break;
			}
			
		}
	}
}

void drawQuad(uint8_t* frameBuffer, Quad* face, PlaydateAPI* pd)
{
	// Minimize float ops, hoist invariants, use integer math where possible
	float x1 = scene->points[face->points[0]].x;
	float y1 = scene->points[face->points[0]].y;
	float x2 = scene->points[face->points[1]].x;
	float y2 = scene->points[face->points[1]].y;
	float x3 = scene->points[face->points[2]].x;
	float y3 = scene->points[face->points[2]].y;
	float x4 = scene->points[face->points[3]].x;
	float y4 = scene->points[face->points[3]].y;

	float u1 = scene->uvs[face->points[0]].u;
	float v1 = scene->uvs[face->points[0]].v;
	float u2 = scene->uvs[face->points[1]].u;
	float v2 = scene->uvs[face->points[1]].v;
	float u3 = scene->uvs[face->points[2]].u;
	float v3 = scene->uvs[face->points[2]].v;
	float u4 = scene->uvs[face->points[3]].u;
	float v4 = scene->uvs[face->points[3]].v;

	// Compute bounding box (use ints for loop bounds)
	float minXf = min4(x1, x2, x3, x4);
	float maxXf = max4(x1, x2, x3, x4);
	float minYf = min4(y1, y2, y3, y4);
	float maxYf = max4(y1, y2, y3, y4);

	int minX = (int)max2(minXf, 0);
	int minY = (int)max2(minYf, 0);
	int maxX = (int)min2(maxXf, LCD_COLUMNS - 1);
	int maxY = (int)min2(maxYf, LCD_ROWS - 1);

	// Precompute triangle setup
	float a12 = y1 - y2, b12 = x2 - x1;
	float a23 = y2 - y3, b23 = x3 - x2;
	float a31 = y3 - y1, b31 = x1 - x3;

	float _a13 = -a31;
	float _a34 = y3 - y4, _b34 = x4 - x3;
	float _a41 = y4 - y1, _b41 = x1 - x4;

	float w12_row = orientation(x1, y1, x2, y2, minX + 0.5f, minY + 0.5f);
	float w23_row = orientation(x2, y2, x3, y3, minX + 0.5f, minY + 0.5f);
	float w31_row = orientation(x3, y3, x1, y1, minX + 0.5f, minY + 0.5f);

	float _w34_row = orientation(x3, y3, x4, y4, minX + 0.5f, minY + 0.5f);
	float _w41_row = orientation(x4, y4, x1, y1, minX + 0.5f, minY + 0.5f);

	float invw123 = 1.f / (w12_row + w23_row + w31_row);
	float invw134 = 1.f / (-w31_row + _w34_row + _w41_row);

	// Precompute UV deltas for scanline stepping
	float deltaU132 = (a23 * u1 + a31 * u2 + a12 * u3) * invw123;
	float deltaV132 = (a23 * v1 + a31 * v2 + a12 * v3) * invw123;
	float deltaU134 = (_a13 * u4 + _a34 * u1 + _a41 * u3) * invw134;
	float deltaV134 = (_a13 * v4 + _a34 * v1 + _a41 * v3) * invw134;

	for(int y = minY; y <= maxY; ++y,
			w23_row += b23, w31_row += b31, w12_row += b12,
			_w34_row += _b34, _w41_row += _b41)
	{
		float w23 = w23_row;
		float w31 = w31_row;
		float w12 = w12_row;
		float _w34 = _w34_row;
		float _w41 = _w41_row;

		uint8_t* target = &frameBuffer[y * 400 + minX];

		float u132 = (w23 * u1 + w31 * u2 + w12 * u3) * invw123;
		float v132 = (w23 * v1 + w31 * v2 + w12 * v3) * invw123;

		float u134 = (-w31 * u4 + _w34 * u1 + _w41 * u3) * invw134;
		float v134 = (-w31 * v4 + _w34 * v1 + _w41 * v3) * invw134;

		int done = 0;
		for(int x = minX; x <= maxX; ++x,
				w23 += a23, w31 += a31, w12 += a12,
				_w34 += _a34, _w41 += _a41, ++target,
				u132 += deltaU132, v132 += deltaV132,
				u134 += deltaU134, v134 += deltaV134)
		{
			if(w31 >= 0.f)
			{
				if(w23 >= 0.f && w12 >= 0.f)
				{
					int _col = (int)u132;
					int _row = (int)v132;
					// Clamp UVs to texture size
					//if(_col >= 0 && _col < 128 && _row >= 0 && _row < 128)
						*target = textureData[_col + _row * 128];
					done = 1;
				}
				else if(done)
				{
					break;
				}
			}
			else //w31 < 0.f
			{
				if(_w34 >= 0.f && _w41 >= 0.f)
				{
					int _col = (int)u134;
					int _row = (int)v134;
					//if(_col >= 0 && _col < 128 && _row >= 0 && _row < 128)
						*target = textureData[_col + _row * 128];
					done = 1;
				}
				else if(done)
				{
					break;
				}
			}
		}
	}
}


void draw(uint8_t* frameBuffer, PlaydateAPI* pd)
{
	uint32_t* fb32 = __frameBuffer__;
	memset(__frameBuffer__, 0, 400 * 240 * sizeof(uint8_t));

	pd->graphics->clear(kColorBlack);
	//pd->graphics->setDrawOffset(0, 0);

	// Draw the current FPS on the screen
	pd->system->drawFPS(0, 0);

	// Draw the scene

	float angle = pd->system->getCrankAngle()/180.f*M_PI;

  if(scene)
  {
		for(size_t i = 0; i < scene->pointCount+1; i++)
		{
			const float factor = 50.f;
			const float dX = 200.f;
			const float dY = 120.f;
			scene->points[i].x = xFactor * (scene->p3d[i].x * cos(angle) - scene->p3d[i].y * sin(angle));
			scene->points[i].x = factor * scene->points[i].x + dX;
			scene->points[i].y = yFactor * (scene->p3d[i].x * sin(angle) + scene->p3d[i].y * cos(angle)) + scene->p3d[i].z;
			scene->points[i].y = factor * scene->points[i].y + dY;
		}

		/*scene->points[10].x = 50.f; scene->points[10].y = 50.f;
		scene->points[11].x = 320.f; scene->points[11].y = 60.f;
		scene->points[12].x = 350.f; scene->points[12].y = 150.f;
		scene->points[13].x = 60.f; scene->points[13].y = 120.f;*/

    for(size_t i = 0; i < scene->faceCount; i++)
    {
      drawFace(__frameBuffer__, &scene->faces[i], pd);
    }

		for(size_t i = 0; i < scene->quadCount; i++)
    {
      drawQuad(__frameBuffer__, &scene->quads[i], pd);
    }

		//drawQuad(__frameBuffer__, &quad, pd);
  }

	uint8_t* target = frameBuffer;
	const uint32_t* source = __frameBuffer__;
	const uint32_t* source2;
	for(size_t row = 0; row < 240; )
	{
		target = frameBuffer + row * LCD_STRIDE;
		source = fb32 + row * 400/4;
		for(size_t col = 0; col < 400/8; ++col, ++target, source += 2)
		{
			source2 = source + 1;
			const uint32_t s1 = *source;
			const uint32_t s2 = *source2;
			*target = ((s1 >> 24) & 1) << 4 |
           ((s1 >> 17) & 1) << 5 |
           ((s1 >> 10)  & 1) << 6 |
           ((s1 >> 3)  & 1) << 7 |
           ((s2 >> 24) & 1) << 0 |
           ((s2 >> 17) & 1) << 1 |
           ((s2 >> 10)  & 1) << 2 |
           ((s2 >> 3)  & 1) << 3;
		}
		++row;
		target = frameBuffer + row * LCD_STRIDE;
		source = fb32 + row * 400/4;
		for(size_t col = 0; col < 400/8; ++col, ++target, source += 2)
		{
			source2 = source + 1;
			const uint32_t s1 = *source;
			const uint32_t s2 = *source2;
			*target = ((s1 >> 27) & 1) << 4 |
           ((s1 >> 16) & 1) << 5 |
           ((s1 >> 9)  & 1) << 6 |
           ((s1 >> 2)  & 1) << 7 |
           ((s2 >> 27) & 1) << 0 |
           ((s2 >> 16) & 1) << 1 |
           ((s2 >> 9)  & 1) << 2 |
           ((s2 >> 2)  & 1) << 3;
		}
		++row;
		target = frameBuffer + row * LCD_STRIDE;
		source = fb32 + row * 400/4;
		for(size_t col = 0; col < 400/8; ++col, ++target, source += 2)
		{
			source2 = source + 1;
			const uint32_t s1 = *source;
			const uint32_t s2 = *source2;
			*target = ((s1 >> 26) & 1) << 4 |
           ((s1 >> 19) & 1) << 5 |
           ((s1 >> 8)  & 1) << 6 |
           ((s1 >> 1)  & 1) << 7 |
           ((s2 >> 26) & 1) << 0 |
           ((s2 >> 19) & 1) << 1 |
           ((s2 >> 8)  & 1) << 2 |
           ((s2 >> 1)  & 1) << 3;
		}
		++row;
		target = frameBuffer + row * LCD_STRIDE;
		source = fb32 + row * 400/4;
		for(size_t col = 0; col < 400/8; ++col, ++target, source += 2)
		{
			source2 = source + 1;
			const uint32_t s1 = *source;
			const uint32_t s2 = *source2;
			*target = ((s1 >> 25) & 1) << 4 |
           ((s1 >> 18) & 1) << 5 |
           ((s1 >> 11)  & 1) << 6 |
           ((s1 >> 0)  & 1) << 7 |
           ((s2 >> 25) & 1) << 0 |
           ((s2 >> 18) & 1) << 1 |
           ((s2 >> 11)  & 1) << 2 |
           ((s2 >> 0)  & 1) << 3;
		}
		++row;
	}
	//pd->system->logToConsole("plop %u", ((uint8_t*)source2 - __frameBuffer__));
}
