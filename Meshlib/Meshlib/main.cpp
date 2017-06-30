#include "simplification.h"
#include "MyMesh.h"


#define BOUNDARY_COST 1.0

using namespace MeshLib;
using namespace std;
CMyMesh mesh;
Simplification simplification;

GLfloat startx, starty;
GLfloat model_angle1 = 0.0, model_angle2 = 0.0, scale = 1.0, eye[3] = { 0.0, 0.0, 10.0 };
GLfloat window_width = 800, window_height = 800;

int Key;
int toggle = 0;
bool left_click = 0, right_click = 0;

bool doEdgeCollapse = false, doVertexSplit = false, doLOD = false;
int  step = 0;



/* window width and height */
int win_width, win_height;
int gButton;
int shadeFlag = 0;

/* rotation quaternion and translation vector for the object */
CQrot       ObjRot(0, 0, 1, 0);
CPoint      ObjTrans(0, 0, 0);

/* arcball object */
CArcball arcball;

int textureFlag = 2;
/* texture id and image */
GLuint texName;
//RgbImage image;


/*! setup the object, transform from the world to the object coordinate system */
void setupObject(void)
{
    double rot[16];
    
    glTranslated(ObjTrans[0], ObjTrans[1], ObjTrans[2]);
    ObjRot.convert(rot);
    glMultMatrixd((GLdouble *)rot);
}

/*! the eye is always fixed at world z = +5 */
void setupEye(void) {
    glLoadIdentity();
    gluLookAt(0, 0, 5, 0, 0, 0, 0, 1, 0);
}

/*! setup light */
void setupLight()
{
    GLfloat lightOnePosition[4] = { 0, 0, 1, 0 };
    glLightfv(GL_LIGHT1, GL_POSITION, lightOnePosition);
}





/*! display call back function
 */
void display()
{
    /* clear frame buffer */
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    setupLight();
    /* transform from the eye coordinate system to the world system */
    //setupEye();
    glPushMatrix();
    /* transform from the world to the ojbect coordinate system */
    setupObject();
    
    glLoadIdentity();
    glScalef(scale, scale, scale);
    gluLookAt(eye[0], eye[1], 10.0,
              eye[0], eye[1],  0.0,
              1.0,    0.0,  0.0 );
    glRotatef(model_angle1, 0, 1, 0);
    glRotatef(model_angle2, 1, 0, 0);
    
    if(doEdgeCollapse){
        simplification.EdgeCollapse();
        doEdgeCollapse = false;
    }
    
    if(doVertexSplit){
        simplification.VertexSplit();
        doVertexSplit = false;
    }
    
    if(doLOD){
        simplification.ControlLevelOfDetail(step);
        doLOD = false;
    }
    
    mesh.Display(toggle);
    
    glPopMatrix();
    glutSwapBuffers();
}

/*! Called when a "resize" event is received by the window. */
void reshape(int w, int h)
{
    float ar;
    //std::cout << "w:" << w << "\th:" << h << std::endl;
    win_width = w;
    win_height = h;
    
    ar = (float)(w) / h;
    glViewport(0, 0, w, h);               /* Set Viewport */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    // magic imageing commands
    gluPerspective(40.0, /* field of view in degrees */
                   ar, /* aspect ratio */
                   1.0, /* Z near */
                   100.0 /* Z far */);
    
    glMatrixMode(GL_MODELVIEW);
    
    glutPostRedisplay();
}

/*! helper function to remind the user about commands, hot keys */
void help()
{
    printf("w  -  Wireframe Display\n");
    printf("+  -  Increase Faces \n");
    printf("-  -  Decrease Faces\n");
    printf("c  -  Edge Collapse\n");
    printf("s  -  Vertex Split\n");
    printf("q  -  quit\n");
}

/*! Keyboard call back function */
void keyBoard(unsigned char key, int x, int y)
{
    switch(key){
        case 'q':
            exit(0);
        case 'c':
            doEdgeCollapse = true;  break;
        case 's':
            doVertexSplit  = true;  break;
        case '-':
            if(step < 200){
                step++;
                doLOD = true;
            }
            break;
        case '+':
            if(step > 0){
                step--;
                doLOD = true;
            }
            break;
        case 'w':
            if(toggle) toggle = 0;
            else       toggle = 1;
            break;
        case '?':
            help();
            break;
    }
    
    glutPostRedisplay();
}


/*! mouse click call back function */
void  mouseClick(int button, int state, int x, int y) {
    /* set up an arcball around the Eye's center
     switch y coordinates to right handed system  */
    
   
        if ( button == GLUT_LEFT_BUTTON ) {
            if ( state == GLUT_DOWN ) {
                left_click = true;
                startx   = x;
                starty   = y;
            }else if (state == GLUT_UP) {
                left_click = false;
            }
        }else{ // button == GLUT_RIGHT_BUTTON
            if ( state == GLUT_DOWN ) {
                right_click = true;
                startx   = x;
                starty   = y;
            }else if (state == GLUT_UP) {
                right_click = false;
            }
        }
    
    
}

/*! mouse motion call back function */
void mouseMove(int x, int y)
{
    if ( left_click && !right_click ) {       // rotation
        model_angle1 += (x - startx);
        model_angle2 += (y - starty);
    }else if( !left_click && right_click ){   // translating
        eye[0] -= (x - startx) / (window_width *0.25);
        eye[1] += (y - starty) / (window_height*0.25);
    }else{ // if( left_click && right_click ) // scaling
        scale -= (y - starty) * 0.01;
    }
    
    startx = x;
    starty = y;
    
    glutPostRedisplay();
    
}


void GLInitt()
{
    float mat_ambient[]   = {0.3, 0.3, 0.3, 1.0};
    float mat_shininess[] = {40.0};
    float mat_specular[]  = {1.0, 1.0, 1.0, 0.0};
    float mat_diffuse[]   = {1.0, 1.0, 1.0, 1.0};
    
    glEnable(GL_DEPTH_TEST);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    // "-near" and "-far" define clipping planes
    glOrtho(-1.1, 1.1, -1.1, 1.1, -10000, 10000);
    
    //gluPerspective(45, 1.0, 1.0, 100);
    
    
    glMatrixMode(GL_MODELVIEW);
    
    static float light0_ambient[] = {0.1, 0.1, 0.1, 1.0};
    static float light0_specular[] = {0.0, 0.0, 0.0, 0.0};
    static float light0_diffuse[] = {0.8, 0.8, 0.8, 0.0};
    static float light0_position[] = {-1000.0, 500.0, 0.0, 0.0};
    static float light1_diffuse[] = {0.8, 0.8, 0.8, 0.0};
    static float light1_position[] = {0.0, 0.0, 1000.0, 0.0};
    static float light2_diffuse[] = {0.8, 0.8, 0.8, 0.0};
    static float light2_position[] = {1000.0, 500.0, 0.0, 0.0};
    static float light_white_diffuse[] = {1.0, 1.0, 1.0, 0.0};
    static float light_white_position[] = {0.0, 00.0, 500.0, 0.0};
    
    glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
    glEnable(GL_LIGHT0);
    
    glLightfv(GL_LIGHT1, GL_AMBIENT, light0_ambient);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
    glLightfv(GL_LIGHT1, GL_SPECULAR, light0_specular);
    glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
    glEnable(GL_LIGHT1);
    /*
     glLightfv(GL_LIGHT2, GL_AMBIENT, light0_ambient);
     glLightfv(GL_LIGHT2, GL_DIFFUSE, light2_diffuse);
     glLightfv(GL_LIGHT2, GL_SPECULAR, light0_specular);
     glLightfv(GL_LIGHT2, GL_POSITION, light2_position);
     glEnable(GL_LIGHT2);
     */
    glLightfv(GL_LIGHT3, GL_AMBIENT, light0_ambient);
    glLightfv(GL_LIGHT3, GL_DIFFUSE, light_white_diffuse);
    glLightfv(GL_LIGHT3, GL_SPECULAR, light0_specular);
    glLightfv(GL_LIGHT3, GL_POSITION, light_white_position);
    
    glEnable(GL_LIGHTING);
    
    glMaterialfv(GL_FRONT, GL_AMBIENT,   mat_ambient);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    glMaterialfv(GL_FRONT, GL_SPECULAR,  mat_specular);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,   mat_diffuse);
    
    glEnable(GL_NORMALIZE);
    glFrontFace(GL_CCW);    // CW: clock wise if we render teapot. Otherwise, use "GL_CCW"
    // glEnable(GL_CULL_FACE);
    // glCullFace(GL_BACK);
}


void init_openGL(int argc, char * argv[])
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_STENCIL);
    glutInitWindowPosition(20, 20);
    glutInitWindowSize(window_width, window_height);
    glutCreateWindow("Mesh Simplification");
    glutDisplayFunc(display);
    glutMouseFunc(mouseClick);
    glutMotionFunc(mouseMove);
    glutKeyboardFunc(keyBoard);
    GLInitt();
    glutMainLoop();
}

/*! main function for viewer
 */
int main(int argc, char * argv[])
{
    
    std::cout<<"Input file:"<<argv[1]<<std::endl;
    
    if (argc != 2|| mesh.ConstructMeshDataStructure(argv[1]) == false)
    {
        std::cout << "Usage: input.off" << std::endl;
        return -1;
    }
    
    simplification.InitSimplification(&mesh);
    
    init_openGL(argc, argv);
    return 0;
}
