// October 28, 2022
// by SM

import std;

#include <GL/glut.h>

import Point3D;
import BsplineSurface;
import RevolvedSurface;

struct graphics
{
    static const double LENGTH;
    static const double oneOverSquareRoot2;
    double tractionAngle{};
    bool leftMouseButtonPressed{};
    bool rightMouseButtonPressed{};
    int lastX{}, lastY{};
    double zoomScale{ LENGTH };
    double aspectRatio{ 1.0 };
    double xTranslation{}, yTranslation{};
    double sNear{ -LENGTH }, sFar{ LENGTH };
    std::array<double, 3> currentVec{}, prevVec{}, rotationAxis{};
    GLdouble mxTransform[4][4]{ {-0.7071, -0.5, 0.5, 0.0}, {0.7071, -0.5, 0.5, 0.0}, {0.0, 0.7071, 0.7071, 0.0}, {0.0, 0.0, 0.0, 1.0} }; // isometric view
};

const double graphics::LENGTH{ 150.0 };
const double graphics::oneOverSquareRoot2{ 1.0 / std::sqrt(2.0) };

graphics g{};

BsplineSurface torus{ 2, 2 }; // first one should be 2; revolution is a circle => degree of two
BsplineSurface gs{ 2, 3 };

void ptTo3DVec(int x, int y, std::array<double, 3>& vec)
{
    // x^2 + y^2 + z^2 == r^2, (r == 1)

    int w{ glutGet(GLUT_WINDOW_WIDTH) };
    int h{ glutGet(GLUT_WINDOW_HEIGHT) };
    //std::cout << std::format("width: {}, heigh: {}\n", w, h);

    vec[0] = 2.0 * x / w - 1.0;
    vec[1] = -2.0 * y / h + 1.0;
    double hypotenuse{ std::hypot(vec[0], vec[1]) };

    if (hypotenuse <= g.oneOverSquareRoot2) // x^2 + y^2 <= r^2 / 2
    {
        vec[2] = std::sqrt(1.0 - hypotenuse * hypotenuse); // z == sqrt(r^2 - (x^2 + y^2))
    }
    else
    {
        vec[2] = 0.5 / hypotenuse; // z == (r^2 / 2) / sqrt(x^2 + y^2)
    }

    hypotenuse = std::hypot(vec[0], vec[1], vec[2]);

    vec[0] /= hypotenuse;
    vec[1] /= hypotenuse;
    vec[2] /= hypotenuse;

    //std::cout << std::format("vector: {}, {}, {}\n", vec[0], vec[1], vec[2]);
}

void onKeyStroke(unsigned char key, int x, int y)
{
    if (key == 'r' || key == 'R')
    {
        g.zoomScale = graphics::LENGTH;
        g.mxTransform[0][0] = -0.7071;
        g.mxTransform[0][1] = -0.5;
        g.mxTransform[0][2] = 0.5;
        g.mxTransform[0][3] = 0.0;
        g.mxTransform[1][0] = 0.7071;
        g.mxTransform[1][1] = -0.5;
        g.mxTransform[1][2] = 0.5;
        g.mxTransform[1][3] = 0.0;
        g.mxTransform[2][0] = 0.0;
        g.mxTransform[2][1] = 0.7071;
        g.mxTransform[2][2] = 0.7071;
        g.mxTransform[2][3] = 0.0;
        g.mxTransform[3][0] = 0.0;
        g.mxTransform[3][1] = 0.0;
        g.mxTransform[3][2] = 0.0;
        g.mxTransform[3][3] = 1.0;
        g.xTranslation = g.yTranslation = 0.0;
        glutPostRedisplay();
    }
}

void onMouseButton(int button, int state, int x, int y)
{
    //std::cout << std::format("button: {}, state: {}, x: {}, y: {}\n", button, state, x, y);
    if (button == 0 && state == 0) // left mouse button pressed
    {
        g.leftMouseButtonPressed = true;
        ptTo3DVec(x, y, g.prevVec);
    }
    else if (button == 0 && state == 1) // left mouse button released
    {
        g.leftMouseButtonPressed = false;
    }
    else if (button == 2 && state == 0) // right mouse button pressed
    {
        g.rightMouseButtonPressed = true;
        g.lastX = x;
        g.lastY = y;
    }
    else if (button == 2 && state == 1) // right mouse button released
    {
        g.rightMouseButtonPressed = false;
    }
    else if (button == 3 && state == 0) // scroll forward
    {
        //std::cout << "scroll forward\n";
        g.zoomScale *= 0.9;
        glutPostRedisplay();
    }
    else if (button == 4 && state == 0) // scroll backward
    {
        //std::cout << "scroll backward\n";
        g.zoomScale *= 1.1;
        glutPostRedisplay();
    }
}

void onMouseDrag(int x, int y)
{
    if (g.leftMouseButtonPressed)
    {
        ptTo3DVec(x, y, g.currentVec);
        //std::cout << std::format("x: {}, y: {}\n", x, y);

        double innerProduct{ g.currentVec[0] * g.prevVec[0] + g.currentVec[1] * g.prevVec[1] + g.currentVec[2] * g.prevVec[2] };
        innerProduct = std::min(innerProduct, 1.0);
        g.tractionAngle = 180.0 * std::acos(innerProduct) / std::numbers::pi; // in degree
        //std::cout << std::format("angle: {}\n", tractionAngle);

        g.rotationAxis[0] = g.prevVec[1] * g.currentVec[2] - g.prevVec[2] * g.currentVec[1];
        g.rotationAxis[1] = g.prevVec[2] * g.currentVec[0] - g.prevVec[0] * g.currentVec[2];
        g.rotationAxis[2] = g.prevVec[0] * g.currentVec[1] - g.prevVec[1] * g.currentVec[0];

        //std::cout << std::format("axis: {}, {}, {}\n", rotationAxis[0], rotationAxis[1], rotationAxis[2]);

        g.prevVec = g.currentVec;

        glutPostRedisplay();
    }
    else if (g.rightMouseButtonPressed)
    {
        g.xTranslation += static_cast<double>(g.lastX - x) * 2.0 * g.zoomScale / glutGet(GLUT_WINDOW_WIDTH);
        g.yTranslation += static_cast<double>(y - g.lastY) * 2.0 * g.zoomScale / glutGet(GLUT_WINDOW_HEIGHT);
        //std::cout << std::format("xTranlation: {}, yTranslation: {}\n", xTranslation, yTranslation);
        g.lastX = x;
        g.lastY = y;
        glutPostRedisplay();
    }
}

void reshape(int x, int y)
{
    g.aspectRatio = static_cast<double>(y) / x; // the inverse of aspect ratio
    glViewport(0, 0, x, y);
}

void display()
{
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-g.zoomScale + g.xTranslation, g.zoomScale + g.xTranslation, -g.zoomScale * g.aspectRatio + g.yTranslation, g.zoomScale * g.aspectRatio + g.yTranslation, g.sNear, g.sFar);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    if (g.leftMouseButtonPressed)
    {
        glPushMatrix();

        glLoadIdentity();
        glRotated(g.tractionAngle, g.rotationAxis[0], g.rotationAxis[1], g.rotationAxis[2]);
        glMultMatrixd(&g.mxTransform[0][0]);
        glGetDoublev(GL_MODELVIEW_MATRIX, &g.mxTransform[0][0]);

        glPopMatrix();
    }

    glMultMatrixd(&g.mxTransform[0][0]);

    glBegin(GL_LINES);
    // x-axis
    glColor3d(1.0, 0.0, 0.0);
    glVertex3d(0.0, 0.0, 0.0);
    glVertex3d(100.0, 0.0, 0.0);

    // y-axis
    glColor3d(0.0, 1.0, 0.0);
    glVertex3d(0.0, 0.0, 0.0);
    glVertex3d(0.0, 100.0, 0.0);

    // z-axis
    glColor3d(0.0, 0.0, 1.0);
    glVertex3d(0.0, 0.0, 0.0);
    glVertex3d(0.0, 0.0, 100.0);
    glEnd();

    // 'x'
    double tx{ 100.0 * g.mxTransform[0][0] + g.mxTransform[3][0] };
    double ty{ 100.0 * g.mxTransform[0][1] + g.mxTransform[3][1] };
    double tz{ 100.0 * g.mxTransform[0][2] + g.mxTransform[3][2] };

    glColor3d(1.0, 0.0, 0.0);

    glPushMatrix();
    glLoadIdentity();
    glTranslated(tx, ty, tz);
    glScaled(0.1, 0.1, 0.1);
    glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, 'x');
    glPopMatrix();

    // 'y'
    tx = 100.0 * g.mxTransform[1][0] + g.mxTransform[3][0];
    ty = 100.0 * g.mxTransform[1][1] + g.mxTransform[3][1];
    tz = 100.0 * g.mxTransform[1][2] + g.mxTransform[3][2];

    glColor3d(0.0, 1.0, 0.0);

    glPushMatrix();
    glLoadIdentity();
    glTranslated(tx, ty, tz);
    glScaled(0.1, 0.1, 0.1);
    glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, 'y');
    glPopMatrix();

    // 'z'
    tx = 100.0 * g.mxTransform[2][0] + g.mxTransform[3][0];
    ty = 100.0 * g.mxTransform[2][1] + g.mxTransform[3][1];
    tz = 100.0 * g.mxTransform[2][2] + g.mxTransform[3][2];

    glColor3d(0.0, 0.0, 1.0);

    glPushMatrix();
    glLoadIdentity();
    glTranslated(tx, ty, tz);
    glScaled(0.1, 0.1, 0.1);
    glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, 'z');
    glPopMatrix();

    glColor3d(1.0, 1.0, 1.0);
    Point3D pt{};

    // draw torus
    for (int u{}; u <= 20; ++u)
    {
        glBegin(GL_LINE_STRIP); // glBegin(GL_POINTS);
        for (int v{}; v <= 20; ++v)
        {
            torus.surfacePointW(u / 20.0, v / 20.0, pt);
            //std::cout << std::format("u: {:15.5f}, v: {:15.5f}, ({:15.5f}, {:15.5f}, {:15.5f})\n", u, v, pt.x, pt.y, pt.z);
            glVertex3d(pt.x, pt.y, pt.z);
        }
        glEnd();
    }

    for (int u{}; u <= 20; ++u)
    {
        glBegin(GL_LINE_STRIP); // glBegin(GL_POINTS);
        for (int v{}; v <= 20; ++v)
        {
            torus.surfacePointW(v / 20.0, u / 20.0, pt);
            glVertex3d(pt.x, pt.y, pt.z);
        }
        glEnd();
    }

    // draw general surface
    for (int u{}; u <= 20; ++u)
    {
        glBegin(GL_LINE_STRIP); // glBegin(GL_POINTS);
        for (int v{}; v <= 20; ++v)
        {
            gs.surfacePointW(u / 20.0, v / 20.0, pt);
            glVertex3d(pt.x, pt.y, pt.z);
        }
        glEnd();
    }

    for (int u{}; u <= 20; ++u)
    {
        glBegin(GL_LINE_STRIP); // glBegin(GL_POINTS);
        for (int v{}; v <= 20; ++v)
        {
            gs.surfacePointW(v / 20.0, u / 20.0, pt);
            glVertex3d(pt.x, pt.y, pt.z);
        }
        glEnd();
    }

    glutSwapBuffers();
}

void makeTorus()
{
    // make torus
    Point3D startPt{ 0,0,0 };
    Point3D directionVec{ 0,0,1 };
    double angle{ 360.0 };
    int m{ 8 };
    std::vector<Point3D> controlPts{ {10, 0, 20}, {10, 0, 30}, {20, 0, 30}, {30, 0, 30}, {30, 0, 20}, {30, 0, 10}, {20, 0, 10}, {10, 0, 10}, {10, 0, 20} };
    std::vector<double> wj{ 1, 0.7071, 1, 0.7071, 1, 0.7071, 1, 0.7071, 1 };
    std::vector<double> V{ 0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1 };

    int n{};
    std::vector<double> U;
    std::vector<std::vector<Point3D>> Pij;
    std::vector<std::vector<double>> wij;

    MakeRevolvedSurface(startPt, directionVec, angle, m, controlPts, wj, n, U, Pij, wij);

    torus.assignControlPoints(Pij);
    torus.assignWeight(wij);

    torus.assignUKnots(U);
    torus.assignVKnots(V);
}

void makeGeneralSurface()
{
    // make torus
    Point3D startPt{ 0, 50, 0 };
    Point3D directionVec{ 1, 0, 0 };
    double angle{ 360.0 };
    int m{ 11 };
    std::vector<Point3D> controlPts{ {20, 40, 0}, {25, 35, 0}, {30, 30, 0}, {40, 30, 0}, {45, 35, 0}, {50, 40, 0}, {55, 45, 0}, {60, 45, 0}, {70, 40, 0}, {75, 35, 0}, {80, 35, 0}, {85, 30, 0} };
    std::vector<double> wj{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    std::vector<double> V{ 0, 0, 0, 0, 0.111111, 0.222222, 0.333333, 0.444444, 0.555555, 0.666666, 0.777777, 0.888888, 1, 1, 1, 1 };
    //std::vector<double> V{ 0, 0, 0, 0.111111, 0.222222, 0.333333, 0.444444, 0.555555, 0.666666, 0.777777, 0.888888, 1, 1, 1 };

    int n{};
    std::vector<double> U;
    std::vector<std::vector<Point3D>> Pij;
    std::vector<std::vector<double>> wij;

    MakeRevolvedSurface(startPt, directionVec, angle, m, controlPts, wj, n, U, Pij, wij);

    gs.assignControlPoints(Pij);
    gs.assignWeight(wij);

    gs.assignUKnots(U);
    gs.assignVKnots(V);
}

int main(int argc, char** argv)
{
    try {
        makeTorus();
        makeGeneralSurface();

        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
        glutInitWindowSize(600, 600);
        glutCreateWindow("Surface of Revolution");
        glutDisplayFunc(display);
        glutReshapeFunc(reshape);
        glutMouseFunc(onMouseButton);
        glutMotionFunc(onMouseDrag);
        glutKeyboardFunc(onKeyStroke);
        glutMainLoop();
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
    catch (...)
    {
        std::cerr << "error\n";
    }
}
