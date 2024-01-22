#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
using namespace std;

#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <GL/glut.h>
#define WIDTH   1024
#define HEIGHT  1024
const int g_TFa = 25;
const int g_TFb = 255;
GLubyte MyTexture[HEIGHT][WIDTH][3]; // 3: rgb
struct voxel {
    float r, g, b, a;
};

voxel volume[256][256][256];
float vol[256][256][256]; // [z][y][x]
float AlphaTable[256];
float ColorTable[256][3];
float sumTable[801];

float* SumTable = nullptr;

float acolor[3] = { 0.0,0.5,0.5 };
float bcolor[3] = { 0.5,0.4,0.6 };
float ccolor[3] = { 0.9,0.0,0.0 };
float dcolor[3] = { 0.0,0.9,0.0 };
void MakeAlphaTable(int a, int b) {
    for (int i = 0; i < 256; i++) {

        if (i < a) {
            AlphaTable[i] = 0;
        }
        else if (i >= a && i < b) {
            AlphaTable[i] = (float)(i - a) / (b - a);
        }
        else if (i >= b) {
            AlphaTable[i] = 1;
        }

        // 원래 함수 분리. 하지만 여기서 한다.
        SumTable = &(sumTable[1]); // 음수 -1 index 사용 가능
        SumTable[-1] = 0.0f;
        for (int j = i; j <= i; j++) {
            SumTable[j] = SumTable[j - 1] + AlphaTable[j];
        }
    }
}
void MakeColorTable(int a, int b, float acolor[3], float bcolor[3], float ccolor[3], float dcolor[3]) {
    for (int i = 0; i < 256; i++) {
        for (int j = 0; j < 3; j++) {

            // 현재 강도에 해당하는 보클 값 결정
            //int voxelValue = static_cast<int>((i - a) / (float)(b - a) * 255.0);
            //// 값을 [0, 255] 범위로 클램핑
            //if (voxelValue < 0) voxelValue = 0;
            //if (voxelValue > 255) voxelValue = 255;

            //// 보클 데이터에서 RGB 값을 사용
            //ColorTable[i][0] = volume[voxelValue][0][0].r ;  // Red
            //ColorTable[i][1] = volume[voxelValue][0][0].g ;  // Green
            //ColorTable[i][2] = volume[voxelValue][0][0].b ;  // Blue
            //// 디버깅을 위한 출력
            //std::cout << "ColorTable[" << i << "]: "
            //   << ColorTable[i][0] << ", "
            //   << ColorTable[i][1] << ", "
            //   << ColorTable[i][2] << ", voxelValue: " << voxelValue << std::endl;
            float t = (float)(i - a) / (b - a);

            if (i < a) {
                ColorTable[i][j] = (1 - t) * acolor[j] + t * bcolor[j];
            }
            else if (i >= a && i < b) {
                ColorTable[i][j] = (1 - t) * bcolor[j] + t * ccolor[j];


            }
            else if (i >= b) {
                ColorTable[i][j] = (1 - t) * ccolor[j] + t * dcolor[j];


            }
        }
    }
    //float ColorTable[256]; // 0~1\
 }
}

bool IsTransparent(int m, int M) { // m==0 인 경우 어쩔?
   //if (m == 0)
   //   return sumTable[M] == 0;
   //else
    return SumTable[M] == SumTable[m - 1]; // -1도 안전함
    for (int i = m; i <= M; i++) {//sumtable 에서 m~M 사이 값변하는지 비교
        if (AlphaTable[i] != 0)
            return false; // 불투명
    }
    return true;
}

class vec {
public:

    float m[3];

    vec() {
        for (int i = 0; i < 3; i++) {
            m[i] = 0;
        }
    }

    vec(float a, float b, float c) {

        m[0] = a;
        m[1] = b;
        m[2] = c;


    }
    vec vec_sub(vec a) {
        vec result;
        for (int i = 0; i < 3; i++) {
            result.m[i] = m[i] - a.m[i];

        }
        return result;
    }

    void vec_print() {
        for (int i = 0; i < 3; i++) {
            printf("(%f )", m[i]);

        }
        printf("\n");

    }

    vec vec_norm() {//벡터 정규화

        float len = sqrtf(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]);

        //예외 처리
        if (len == 0)
            len = 1;
        for (int i = 0; i < 3; i++) {

            m[i] = m[i] / len;
        }



        return *this;

    }
    vec operator+(const vec& other) const {
        return vec(m[0] + other.m[0], m[1] + other.m[1], m[2] + other.m[2]);
    }

    vec operator*(float scalar) const {
        return vec(m[0] * scalar, m[1] * scalar, m[2] * scalar);
    }
    vec operator -(vec v) {
        return vec(m[0] - v.m[0], m[1] - v.m[1], m[2] - v.m[2]);
    }
    vec operator +(vec v) {
        return vec(m[0] + v.m[0], m[1] + v.m[1], m[2] + v.m[2]);
    }

    float operator*(vec ex) const {
        return m[0] * ex.m[0] + m[1] * ex.m[1] + m[2] * ex.m[2];
    }

    float dot(vec x) {
        return m[0] * x.m[0] + m[1] * x.m[1] + m[2] * x.m[2];
    }





    vec vec_cross(vec a) { // c = a * b  외적
       //cross product가 안나오는 경우 -> 두 벡터가 붙어있는 경우 답이 무한히 많아짐
        vec result;
        result.m[0] = m[1] * a.m[2] - m[2] * a.m[1];
        result.m[1] = m[2] * a.m[0] - m[0] * a.m[2];
        result.m[2] = m[0] * a.m[1] - m[1] * a.m[0];

        return result;
    }

};
//float eye[3] = { 0,160,128 };
//float at[3] = { 128,128,112.5 };
//float dir[3];
//float up[3] = { 0,-1,0 };
//
//float u[3], v[3], w[3]; // 카메라 기준의  x,y,z축 방향
vec eye(0, 200, 256);
vec at(128, 128, 128);
vec dir;
vec up(0, -1, 0);
vec w;
vec u;
vec v;
vec rayStart;
void vec_print(float x[3]) {
    printf("(%f, %f, %f)\n", x[0], x[1], x[2]);
}
// (0,0,0)~from 까지의 벡터를 normalize 해서 to 에 저장
void vec_norm(float from[3], float to[3]) {
    float len = sqrtf(from[0] * from[0] + from[1] * from[1] + from[2] * from[2]);
    if (len == 0)
        len = 1;// 예외적으로 그냥 복사해줌
    to[0] = from[0] / len;
    to[1] = from[1] / len;
    to[2] = from[2] / len;
}
void vec_sub(float c[3], float a[3], float b[3]) {
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    c[2] = a[2] - b[2];
}
void vec_cross(float c[3], float a[3], float b[3]) { // c = a * b
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

void samplingRGB(float x, float y, float z, float* color) {
    if (x >= 255 || y >= 255 || z >= 255 || x < 0 || y < 0 || z < 0) {
        color[0] = color[1] = color[2] = 0;
        return;
    }
    int ix = (int)x;
    float wx = x - ix;
    int iy = (int)y;
    float wy = y - iy;
    int iz = (int)z;
    float wz = z - iz;

    for (int i = 0; i < 3; i++) { // 각 색상 채널에 대해 별도로 보간
        float a = (i == 0 ? volume[iz][iy][ix].r : (i == 1 ? volume[iz][iy][ix].g : volume[iz][iy][ix].b));
        float b = (i == 0 ? volume[iz][iy][ix + 1].r : (i == 1 ? volume[iz][iy][ix + 1].g : volume[iz][iy][ix + 1].b));
        float c = (i == 0 ? volume[iz][iy + 1][ix].r : (i == 1 ? volume[iz][iy + 1][ix].g : volume[iz][iy + 1][ix].b));
        float d = (i == 0 ? volume[iz][iy + 1][ix + 1].r : (i == 1 ? volume[iz][iy + 1][ix + 1].g : volume[iz][iy + 1][ix + 1].b));
        float e = (i == 0 ? volume[iz + 1][iy][ix].r : (i == 1 ? volume[iz + 1][iy][ix].g : volume[iz + 1][iy][ix].b));
        float f = (i == 0 ? volume[iz + 1][iy][ix + 1].r : (i == 1 ? volume[iz + 1][iy][ix + 1].g : volume[iz + 1][iy][ix + 1].b));
        float g = (i == 0 ? volume[iz + 1][iy + 1][ix].r : (i == 1 ? volume[iz + 1][iy + 1][ix].g : volume[iz + 1][iy + 1][ix].b));
        float h = (i == 0 ? volume[iz + 1][iy + 1][ix + 1].r : (i == 1 ? volume[iz + 1][iy + 1][ix + 1].g : volume[iz + 1][iy + 1][ix + 1].b));

        color[i] = a * (1 - wx) * (1 - wy) * (1 - wz) + b * (wx) * (1 - wy) * (1 - wz) +
            c * (1 - wx) * wy * (1 - wz) + d * wx * wy * (1 - wz) + e * (1 - wx) * (1 - wy) * wz
            + f * (wx) * (1 - wy) * wz + g * (1 - wx) * wy * wz + h * wx * wy * wz;
    }
}



float sampling(float x, float y, float z) { // x=3.7 , ix = 3
   // 힌트, 좌표 벗어나서 메모리 참조 오류 해결할 것.
    if (x >= 255 || y >= 255 || z >= 255 ||
        x < 0 || y < 0 || z < 0)
        return 0;
    int ix = (int)x;
    float wx = x - ix;
    int iy = (int)y;
    float wy = y - iy;
    int iz = (int)z;
    float wz = z - iz;


    float a = vol[iz][iy][ix];
    float b = vol[iz][iy][ix + 1];
    float c = vol[iz][iy + 1][ix];
    float d = vol[iz][iy + 1][ix + 1];
    float e = vol[iz + 1][iy][ix];
    float f = vol[iz + 1][iy][ix + 1];
    float g = vol[iz + 1][iy + 1][ix];
    float h = vol[iz + 1][iy + 1][ix + 1];

    float res = a * (1 - wx) * (1 - wy) * (1 - wz) + b * (wx) * (1 - wy) * (1 - wz) +
        c * (1 - wx) * wy * (1 - wz) + d * wx * wy * (1 - wz) + e * (1 - wx) * (1 - wy) * wz
        + f * (wx) * (1 - wy) * wz + g * (1 - wx) * wy * wz + h * wx * wy * wz;


    //if (res > 255) res = 255;
    return res; // 77.285

//return vol[iz][iy][ix];

//float a = vol[iz][iy][ix], float b=[iz][iy][ix+1], float c, float d
//float res = a* (1 - wx)* (1 - wy) + b*(wx) * (1 - wy) +
//   c * (1 - wx) * wy + d * wx * wy;
//return res; // 77.285
}

float maxr = 0;
float maxg = 0;
float maxb = 0;
float maxa = 0;
float maxvol = 0;
void FillMyTexture() {
    // FILE* fp = fopen("256x256x256_0.den", "rb");
    FILE* fp = fopen("rgba_data.den", "rb");
    fread(volume, sizeof(voxel), 256 * 256 * 256, fp);
    fclose(fp);
    for (int iy = 0; iy < 256; iy++) {
        for (int ix = 0; ix < 256; ix++) {
            for (int iz = 0; iz < 256; iz++) {
                vol[iy][ix][iz] = volume[iy][ix][iz].a * 10;
                //추후 코드 수정
                maxr = __max(maxr, volume[iy][ix][iz].r);
                maxg = __max(maxg, volume[iy][ix][iz].g);
                maxb = __max(maxb, volume[iy][ix][iz].b);
                maxa = __max(maxa, volume[iy][ix][iz].a);
                maxvol = __max(maxvol, vol[iy][ix][iz]);

            }
        }
    }
    printf("fun get_rgba_on_grid max r:!! %f !! \n", maxr);
    printf("fun get_rgba_on_grid max g:!! %f !! \n", maxg);
    printf("fun get_rgba_on_grid max b:!! %f !! \n", maxb);
    printf("fun get_rgba_on_grid max a:!! %f !! \n", maxa);
    printf("fun vol max a:!! %f !! \n", maxvol);

    // MPR: MultiPlanar Reformation
    // MIP : Maximum Intensity Projection
    float max = 0;

    for (int iy = 0; iy < HEIGHT; iy++) {
        for (int ix = 0; ix < WIDTH; ix++) {
            int dx = (ix - 512) / 2; // 중심점 기준으로 영상 u축 방향으로 몇칸 움직였나
            int dy = (iy - 512) / 2; // 중심점 기준으로 영상 v축 방향으로 몇칸 움직였나

            //vec3 rs = eye + u * dx + v * dy;
            rayStart = eye + u * dx + v * dy;
            int maxv = 0;
            float alphasum = 0;
            float colorsum[3];
            for (int t = 0; t < 3; t++) {
                colorsum[t] = 0;
            }
            for (int t = 0; t < 255; t++) {
                // vec3 s = rs + w*t;
                vec s;
                s.m[0] = rayStart.m[0] + w.m[0] * t;
                s.m[1] = rayStart.m[1] + w.m[1] * t;
                s.m[2] = rayStart.m[2] + w.m[2] * t;
                //float vz = t, vy = iy, vx = ix; // ity % 256, vy = ix % 256, vx = t;
                //if (vz >= 225)
                //   continue;
                GLubyte Intensity = sampling(s.m[0], s.m[1], s.m[2]); // samplig(s);

                float color[3];
                samplingRGB(s.m[0], s.m[1], s.m[2], color);
                max = __max(max, Intensity);

                float alpha = AlphaTable[Intensity];
                if (alpha == 0) {
                    continue;
                }
                float dx = (sampling(s.m[0] + 1, s.m[1], s.m[2]) - sampling(s.m[0] - 1, s.m[1], s.m[2])) * 0.5;
                float dy = (sampling(s.m[0], s.m[1] + 1, s.m[2]) - sampling(s.m[0], s.m[1] - 1, s.m[2])) * 0.5;
                float dz = (sampling(s.m[0], s.m[1], s.m[2] + 1) - sampling(s.m[0], s.m[1], s.m[2] - 1)) * 0.5;
                vec N;
                N.m[0] = dx;
                N.m[1] = dy;
                N.m[2] = dz;
                N.vec_norm();
                // alpha blending

                //for (int i = 0; i < 3; i++) {
                //    color[i] = ColorTable[Intensity][i]; // ColorTable에는 r,g,b가 들어있다. color[3] = {0.4, 0.2, 0.1}; // 누르스름 주황색
                //}
                float Ia = 0.3, Id = 0.6, Is = 0.3; // 고정, Ia[3] = {0.3 , 0.3, 0.3}; // r,g,b동일 느낌
                float Ka[3], Kd[3];
                for (int i = 0; i < 3; i++) {

                    Ka[i] = color[i];
                    Kd[i] = color[i];

                }
                float Ks = 1.0;
                vec L(w); // 광부 조명 모델(편리)
                L = L * -1.0;

                float NL = N.dot(L);
                if (NL < 0)
                    NL = -NL;
                vec H = eye; // 광부모델의 경우만 가능 원래는 (L + eye).normalize();
                float NH = N.dot(H); // RV
                if (NH < 0)
                    NH = -NH;
                vec R = N * (2 * NL) - L;
                R.vec_norm();
                float RV = fabs(R.dot(w));
                const float po = 10;
                // K는 물체색상, 그러므로 color[3]을 사용.
                float I[3];
                float e[3];
                for (int i = 0; i < 3; i++) {

                    I[i] = Ia * Ka[i] + Id * Kd[i] * NL + Is * Ks * pow(RV, po);
                    e[i] = I[i] * alpha;
                    colorsum[i] = colorsum[i] + e[i] * (1 - alphasum);

                }


                alphasum = alphasum + alpha * (1 - alphasum);
                // 속도향상: 화질을 저하없이 (최소한으로하면서) 연산을 줄이는 방법 연구
                if (alphasum > 0.99) {
                    break;
                }

            }

            for (int i = 0; i < 3; i++) {

                if (colorsum[i] > 1)
                    colorsum[i] = 1;
            }

            //GLubyte Intensity = vol[(int)vz][(int)vy][(int)vx]; // ((s + t) % 2) * 255;    //0는 흑색, 255는 백색
            //maxv = __max(maxv, Intensity);

         //Intensity = rand()%255;
            MyTexture[iy][ix][0] = colorsum[0] * 255.0;         //Red 값에 할당
            MyTexture[iy][ix][1] = colorsum[1] * 255.0;;             //Green 값에 할당
            MyTexture[iy][ix][2] = colorsum[2] * 255.0;
            //Blue 값에 할당

        }
    }
    printf(" max sss : %f\n", max);
}

void MyInit() {
    glClearColor(0.0, 0.0, 0.0, 0.0);
    FillMyTexture();
    glTexImage2D(GL_TEXTURE_2D, 0, 3, WIDTH, HEIGHT, 0, GL_RGB,
        GL_UNSIGNED_BYTE, &MyTexture[0][0][0]);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP); // GL_REPEAT
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP); // GL_REPEAT
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
    glEnable(GL_TEXTURE_2D);
}

void MyDisplay() {
    glClear(GL_COLOR_BUFFER_BIT);
    glBegin(GL_QUADS);
    // (0, 0) ~ (3, 3)
    glTexCoord2f(0.0, 1.0); glVertex3f(-1.0, -1.0, 0.0);
    glTexCoord2f(0.0, 0.0); glVertex3f(-1.0, 1.0, 0.0);
    glTexCoord2f(1.0, 0.0); glVertex3f(1.0, 1.0, 0.0);
    glTexCoord2f(1.0, 1.0); glVertex3f(1.0, -1.0, 0.0);
    glEnd();
    glutSwapBuffers();
}
int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitWindowSize(1024, 1024);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutCreateWindow("openGL Sample Program");

    dir = at.vec_sub(eye);
    //vec_norm(dir, w);//dir을 w로 단위벡터화
    w = dir.vec_norm();

    //up벡터가 주어졌을때, w방향으로 정사영(그림자)를 내려서 만나는 방향이 v[3]가 되도록 계산하라
    //v[3]가 주어진다면 u[3]은 자동적으로 결정
    //쉽게 풀기 위해 외적 활용(vec_cross)
    //오른손 법칙에 의해 방향이 결정
    //up, w를 외적 시키면 u 즉 x축
    //w, u를 외적 시키면 v 즉 y축 => up벡터의 수선의 발
    u = up.vec_cross(w); //xyzxyzxyzxyz  u는 x축
    u.vec_norm();
    v = w.vec_cross(u);//외적을 이용해 y 축 구함
    v.vec_norm();

    printf("--uvw---\n");
    u.vec_print(); // x축
    v.vec_print(); // y축
    w.vec_print(); // z축
    MakeAlphaTable(g_TFa, g_TFb);//투명도
    //MakeColorTable(0, 255, acolor, bcolor, ccolor, dcolor);//색
    float max = 0;
    for (int i = 0; i < 255; i++) {
        for (int j = 0; j < 3; j++) {
            max = __max(max, ColorTable[i][j]);
        }
    }
    printf("color max : %f\n", max);

    MyInit();

    glutDisplayFunc(MyDisplay);
    glutMainLoop();
    return 0;
}