#ifndef _MY_MESH_
#define _MY_MESH_

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <queue>
#include <deque>
#include <stack>
#include <GLUT/glut.h>
#include <math.h>
#include <stdio.h>
#include <cstdio>
#include <stdlib.h>

#include "core/Mesh/Vertex.h"
#include "core/Mesh/Edge.h"
#include "core/Mesh/Face.h"
#include "core/Mesh/HalfEdge.h"
#include "core/Mesh/BaseMesh.h"
#include "core/viewer/Arcball.h"                           /*  Arc Ball  Interface         */
#include "core/Mesh/boundary.h"
#include "core/Mesh/iterators.h"
#include "core/Parser/parser.h"



#ifndef M_PI
#define M_PI 3.141592653589793238
#define EPSILON 1.0e-6
#endif

namespace MeshLib
{
    class CMyVertex;
    class CMyEdge;
    class CMyFace;
    class CMyHalfEdge;
    typedef std::list<CMyVertex>::iterator  VertexIter;
    typedef std::list<CMyFace>::iterator    FaceIter;
    typedef std::list<CMyEdge>::iterator    EdgeIter;
    
    class CMyVertex : public CVertex
    {
    public:
        CMyVertex() : m_rgb(1,1,1) {};
        CMyVertex(double *coord_in, int n){
            for(int i=0; i<3; i++) coord[i] = coord_in[i];
            id = n;
            neighborHe = NULL;
            isBoundary = false;
            isActive = true;
        }
        ~CMyVertex() {};
        
        CMyHalfEdge *neighborHe;
        double coord[3];
        double normal[3];
        double Q[10];
        int id;
        bool isBoundary;
        bool isActive;
        
        void _from_string() ;
        void _to_string();
        CPoint & rgb() { return m_rgb; };
    protected:
        CPoint m_rgb;
    };
    
    inline void CMyVertex::_from_string()
    {
        CParser parser(m_string);
        for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
        {
            CToken * token = *iter;
            if (token->m_key == "uv") //CPoint2
            {
                token->m_value >> m_uv;
            }
            if (token->m_key == "rgb") // CPoint
            {
                token->m_value >> m_rgb;
            }
        }
    }
    
    inline void CMyVertex::_to_string()
    {
        CParser parser(m_string);
        parser._removeToken("uv");
        
        parser._toString(m_string);
        std::stringstream iss;
        
        iss << "uv=(" << m_uv[0] << " " << m_uv[1] << ")";
        
        if (m_string.length() > 0)
        {
            m_string += " ";
        }
        m_string += iss.str();
    }
    
    class CMyEdge : public CEdge
    {
    public:
        CMyHalfEdge *halfedge[2];
        CMyEdge() :m_sharp(false) {};
        CMyEdge(CMyHalfEdge *he0, CMyHalfEdge *he1, int n){
            halfedge[0] = he0;
            halfedge[1] = he1;
            id = n;
            isActive = true;
        };

        ~CMyEdge() {};
        
        void _from_string();
        bool & sharp() { return m_sharp; };
        int id;
        bool m_sharp;
        int ect_id;
        bool isActive;
        
    };
    
    inline void CMyEdge::_from_string()
    {
        CParser parser(m_string);
        for (std::list<CToken*>::iterator iter = parser.tokens().begin(); iter != parser.tokens().end(); ++iter)
        {
            CToken * token = *iter;
            if (token->m_key == "sharp")
            {
                m_sharp = true;
            }
        }
    }
    
    
    class CMyHalfEdge : public CHalfEdge
    {
    public:
        CMyHalfEdge *next, *prev;
        CMyHalfEdge *mate;
        
        VertexIter vertex;
        FaceIter   face;
        EdgeIter   edge;
        
        CMyHalfEdge(){
            next = prev = NULL;
            mate = NULL;
        }
    };
    
    
    class CMyFace : public CFace
    {
    public:
        CMyFace(){}
        CMyFace(VertexIter v0, VertexIter v1, VertexIter v2, int n){
            halfedge[0].vertex = v0;
            halfedge[1].vertex = v1;
            halfedge[2].vertex = v2;
            id = n;
            isActive = true;
        }
        CPoint & normal() { return m_normal; };
        CMyHalfEdge   halfedge[3];
        
        double ormal[3];
        double area;
        int id;
        bool isActive;
    protected:
        CPoint m_normal;
    };
    
    
    
    template<typename V, typename E, typename F, typename H>
    class MyMesh
    {
    public:
        typedef CBoundary<V, E, F, H>					CBoundary;
        typedef CLoop<V, E, F, H>						CLoop;
        
        typedef MeshVertexIterator<V, E, F, H>			MeshVertexIterator;
        typedef MeshEdgeIterator<V, E, F, H>			MeshEdgeIterator;
        typedef MeshFaceIterator<V, E, F, H>			MeshFaceIterator;
        typedef MeshHalfEdgeIterator<V, E, F, H>		MeshHalfEdgeIterator;
        
        typedef VertexVertexIterator<V, E, F, H>		VertexVertexIterator;
        typedef VertexEdgeIterator<V, E, F, H>			VertexEdgeIterator;
        typedef VertexFaceIterator<V, E, F, H>			VertexFaceIterator;
        typedef VertexInHalfedgeIterator<V, E, F, H>	VertexInHalfedgeIterator;
        typedef VertexOutHalfedgeIterator<V, E, F, H>	VertexOutHalfedgeIterator;
        
        typedef FaceVertexIterator<V, E, F, H>			FaceVertexIterator;
        typedef FaceEdgeIterator<V, E, F, H>			FaceEdgeIterator;
        typedef FaceHalfedgeIterator<V, E, F, H>		FaceHalfedgeIterator;
        
        
        std::list<CMyVertex> vertices;
        std::list<CMyFace>   faces;
        std::list<CMyEdge>   edges;
        
        
        int n_vertices, n_faces, n_edges;
        
        
        MyMesh(){
            n_vertices = n_faces = n_edges = 0;
        }
        bool Read(char *filename);
        bool ConstructMeshDataStructure(char *filename);
        void output_mesh_info();
        void test_iterator();
        void AddEdgeInfo();
        void MakeCircularList(FaceIter &fi);
        void AssignFaceNormal(FaceIter &fi);
        void AssignVertexNormal(VertexIter &vi);
        void Display(int mode);
    };
    
    typedef MyMesh<CMyVertex, CMyEdge, CMyFace, CMyHalfEdge> CMyMesh;
    
    template<typename V, typename E, typename F, typename H>
    void MeshLib::MyMesh<V, E, F, H>::output_mesh_info()
    {
        int nv = this->numVertices();
        int ne = this->numEdges();
        int nf = this->numFaces();
        
        std::cout << "#V=" << nv << "  ";
        std::cout << "#E=" << ne << "  ";
        std::cout << "#F=" << nf << "  ";
        
        int euler_char= nv - ne + nf;
        std::cout << "Euler's characteristic=" << euler_char << "  ";
        
        CBoundary boundary(this);
        std::vector<CLoop*> & loops = boundary.loops();
        int nb = loops.size();
        
        int genus = (2 - (euler_char + nb)) / 2;
        std::cout << "genus=" << genus << std::endl;
    }
    
    template<typename V, typename E, typename F, typename H>
    void MyMesh<V, E, F, H>::test_iterator()
    {
        for (MeshVertexIterator viter(this); !viter.end(); ++viter)
        {
            V * pV = *viter;
            // you can do something to the vertex here
            // ...
            
            for (VertexEdgeIterator veiter(pV); !veiter.end(); ++veiter)
            {
                E * pE = *veiter;
                // you can do something to the neighboring edges with CCW
                // ...
            }
            
            for (VertexFaceIterator vfiter(pV); !vfiter.end(); ++vfiter)
            {
                F * pF = *vfiter;
                // you can do something to the neighboring faces with CCW
                // ...
            }
            
            for (VertexInHalfedgeIterator vhiter(this, pV); !vhiter.end(); ++vhiter)
            {
                H * pH = *vhiter;
                // you can do something to the incoming halfedges with CCW
                // ...
            }
        }
        
        for (MeshEdgeIterator eiter(this); !eiter.end(); ++eiter)
        {
            E * pE = *eiter;
            // you can do something to the edge here
            // ...
        }
        
        for (MeshFaceIterator fiter(this); !fiter.end(); ++fiter)
        {
            F * pF = *fiter;
            // you can do something to the face here
            // ...
        }
        
        //there are some other iterators which you can find them in class MyMesh
        
        std::cout << "Iterators test OK.\n";
    }
    
    extern void GLInit();
    extern void CrossProduct( double *a, double *b, double *c );
    extern void Normalize(double *a);
    extern double DotProduct(double *a, double *b);
    extern double DotProduct4D(double *a, double *b);
    extern void GetEdgeVector(double *v1, double *v2, double *edge_vector);
    extern void GetArea(double *normal, double &area);
    extern double GetDistance(double *a, double *b);
    extern double GetLength(double *a);
    extern void Swap(double &a, double &b);
    extern bool SolveLinearSystem(double (*matrix)[4], double *rhs, double *solution);


inline void CrossProduct(double *a, double *b, double *c)
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

inline void Normalize(double *a)
{
    double length = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    
    if(!length) return;
    
    double inv_length = 1.0/length;
    
    a[0] *= inv_length;
    a[1] *= inv_length;
    a[2] *= inv_length;
}

inline double DotProduct(double *a, double *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline double DotProduct4D(double *a, double *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
}

inline void GetArea(double *normal, double &area)
{
    area = 0.5 * sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
}


inline double GetDistance(double *a, double *b)
{
    return sqrt( pow(a[0]-b[0], 2.0) + pow(a[1]-b[1], 2.0) + pow(a[2]-b[2], 2.0) );
}

inline double GetLength(double *a)
{
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}


inline void Swap(double &a, double &b)
{
    double t = a;  a = b;  b = t;
}

inline bool SolveLinearSystem(double (*matrix)[4], double *rhs, double *solution)
{
    for(int i = 0; i <= 3; i++){
        double maxAbsoluteValue = -1.0;
        int    pivot_index;
        
        for(int j = i; j <= 3; j++){
            if(fabs(matrix[j][i]) > maxAbsoluteValue){
                maxAbsoluteValue = fabs(matrix[j][i]);
                pivot_index      = j;
            }
        }
        
        if(maxAbsoluteValue < 1.0e-6) return false;
        
        for(int j = i; j <= 3; j++) Swap(matrix[i][j], matrix[pivot_index][j]);
        Swap(rhs[i], rhs[pivot_index]);
        
        double scale = 1.0 / matrix[i][i];
        
        for(int j = i+1; j <= 3; j++){
            double pivot = -matrix[j][i]*scale;
            
            for(int k = 0;   k <= i; k++) matrix[j][k] = 0.0;
            for(int k = i+1; k <= 3; k++){
                if(fabs(matrix[i][k]) > 1.0e-6)
                    matrix[j][k] += matrix[i][k] * pivot;
                else
                    break;
            }
            
            rhs[j] += rhs[i] * pivot;
        }
    }
    
    for(int i = 3; i >= 0; i--){
        solution[i] = 0.0;
        for(int j = i+1; j <= 3; j++) solution[i] += solution[j]*matrix[i][j];
        
        solution[i] = (rhs[i]-solution[i])/matrix[i][i];
    }
    
    
    return true;
}
    
    template<> void CMyMesh::AssignVertexNormal(VertexIter &vi)
    {
        bool isBoundaryVertex = false;
        
        vi->normal[0] = vi->normal[1] = vi->normal[2] = 0.0;
        double cumulativeArea = 0.0;
        
        CMyHalfEdge *hep = vi->neighborHe;
        do{
            FaceIter fi = hep->face;
            
            vi->normal[0] += fi->ormal[0]*fi->area;
            vi->normal[1] += fi->ormal[1]*fi->area;
            vi->normal[2] += fi->ormal[2]*fi->area;
            cumulativeArea += fi->area;
            
            hep = hep->prev->mate;
            
            if(hep == NULL){
                isBoundaryVertex = true;
                break;
            }
        }while(hep != vi->neighborHe);
        
        if(isBoundaryVertex){
            CMyHalfEdge *hep = vi->neighborHe->mate;
            while(hep != NULL){
                FaceIter fi = hep->face;
                
                vi->normal[0] += fi->ormal[0]*fi->area;
                vi->normal[1] += fi->ormal[1]*fi->area;
                vi->normal[2] += fi->ormal[2]*fi->area;
                cumulativeArea += fi->area;
                
                hep = hep->next->mate;
            }
        }
        
        double invCumulativeArea = 1.0 / cumulativeArea;
        
        vi->normal[0] *= invCumulativeArea;
        vi->normal[1] *= invCumulativeArea;
        vi->normal[2] *= invCumulativeArea;
    }
    
    template<> void CMyMesh::AssignFaceNormal(FaceIter &fi)
    {
        double vec1[3], vec2[3];
        
        for(int i = 0; i < 3; i++){
            vec1[i] = fi->halfedge[1].vertex->coord[i] - fi->halfedge[0].vertex->coord[i];
            vec2[i] = fi->halfedge[2].vertex->coord[i] - fi->halfedge[0].vertex->coord[i];
        }
        
        CrossProduct(vec1, vec2, fi->ormal);
        GetArea(fi->ormal, fi->area);
        Normalize(fi->ormal);
    }
    
    template<> void CMyMesh::MakeCircularList(FaceIter &fi)
    {
        fi->halfedge[0].next = &(fi->halfedge[1]);
        fi->halfedge[1].next = &(fi->halfedge[2]);
        fi->halfedge[2].next = &(fi->halfedge[0]);
        
        fi->halfedge[0].prev = &(fi->halfedge[2]);
        fi->halfedge[1].prev = &(fi->halfedge[0]);
        fi->halfedge[2].prev = &(fi->halfedge[1]);
    }

    template<> void CMyMesh::AddEdgeInfo()
    {
        std::vector< std::vector<FaceIter> > Ring(n_vertices);
        
        for(FaceIter fi = faces.begin(); fi != faces.end(); fi++){
            
            MakeCircularList(fi);
            
            for(int i = 0; i < 3; i++){
                fi->halfedge[i].face = fi;
                fi->halfedge[i].vertex->neighborHe = &(fi->halfedge[i]);
                
                Ring[fi->halfedge[i].vertex->id].push_back(fi);
            }
        }
        
            for(unsigned int i = 0; i < Ring.size(); i++){
            
            for(unsigned int j = 0; j < Ring[i].size(); j++){
                
                CMyHalfEdge   *candidate_he = nullptr;
                VertexIter candidate_vertex;
                
                for(int m = 0; m < 3; m++){
                    if(Ring[i][j]->halfedge[m].vertex->id == i){
                        candidate_he     = &(Ring[i][j]->halfedge[m]);
                        candidate_vertex = Ring[i][j]->halfedge[m].next->vertex;
                        break;
                    }
                }
                
                for(unsigned int k = 0; k < Ring[i].size(); k++){
                    
                    if(j==k) continue;
                    
                    for(int m = 0; m < 3; m++){
                        if( Ring[i][k]->halfedge[m].vertex == candidate_vertex &&
                           Ring[i][k]->halfedge[m].next->vertex->id == i ){
                            candidate_he->mate = &(Ring[i][k]->halfedge[m]);
                            Ring[i][k]->halfedge[m].mate = candidate_he;
                            break;
                        }
                    }
                    
                }
            }
        }
        
        for(FaceIter fi = faces.begin(); fi != faces.end(); fi++){
            
            for(int i = 0; i < 3; i++){
                if( fi->halfedge[i].mate == NULL ||
                   fi->halfedge[i].vertex->id < fi->halfedge[i].mate->vertex->id ){
                    edges.push_back( CMyEdge( &(fi->halfedge[i]), fi->halfedge[i].mate, n_edges++) );
                }
                if( fi->halfedge[i].mate == NULL ) fi->halfedge[i].vertex->isBoundary = true;
            }
        }
        
        for(EdgeIter ei = edges.begin(); ei != edges.end(); ei++){
            ei->halfedge[0]->edge = ei;
            if(ei->halfedge[1] != NULL) ei->halfedge[1]->edge = ei;
        }
        
        std::cerr << "Number of edges     "  << n_edges << std::endl;
        
        for(FaceIter   fi = faces.begin();    fi != faces.end();    fi++) AssignFaceNormal(fi);
        for(VertexIter vi = vertices.begin(); vi != vertices.end(); vi++) AssignVertexNormal(vi);
    }
    
    
    template<> bool CMyMesh::Read(char *filename)
    {
        FILE *fp;
        
        if((fp = fopen(filename, "r")) == NULL ){
            std::cerr << "file cannot be read.\n";
            return false;
        }
        
        char buf[512];
        
        fgets(buf, 512, fp);
        fgets(buf, 512, fp);
        sscanf(buf, "%d%d", &n_vertices, &n_faces);
        
        std::vector<VertexIter> vertex_iterator;
        int v_id = 0, f_id = 0;
        
        for(int i = 0; i < n_vertices; i++){
            double coord_in[3];
            
            fgets(buf, 512, fp);
            sscanf(buf, "%lf%lf%lf", &coord_in[0], &coord_in[1], &coord_in[2]);
            
            vertices.push_back( CMyVertex(coord_in, v_id++) );
            vertex_iterator.push_back( --(vertices.end()) );
        }
        
        for(int i = 0; i < n_faces; i++){
            int v_id[3], dummy;
            
            fgets(buf, 512, fp);
            sscanf(buf, "%d%d%d%d", &dummy, &v_id[0], &v_id[1], &v_id[2]);
            
            faces.push_back( CMyFace(vertex_iterator[ v_id[0] ], vertex_iterator[ v_id[1] ], vertex_iterator[ v_id[2] ], f_id++) );
        }
        
        std::cerr << "Read file done...\n";
        
        fclose(fp);
        
        
        std::cerr << "Number of vertices  " << n_vertices << std::endl;
        std::cerr << "Number of faces     " << n_faces    << std::endl;
        
        
        double range_min[3] = {  1.0e6,  1.0e6,  1.0e6, };
        double range_max[3] = { -1.0e6, -1.0e6, -1.0e6, };
        double center[3];
        
        for(VertexIter vi = vertices.begin(); vi != vertices.end(); vi++){
            for(int i = 0; i < 3; i++){
                if(vi->coord[i] < range_min[i])	range_min[i] = vi->coord[i];
                if(vi->coord[i] > range_max[i])	range_max[i] = vi->coord[i];
            }
        }
        
        for(int i = 0; i < 3; i++) center[i] = (range_min[i] + range_max[i])*0.5;
        
        double largest_range = -1.0;
        
        for(int i = 0; i < 3; i++){
            if(largest_range < range_max[i]-range_min[i]) largest_range = range_max[i]-range_min[i];
        }
        
        double scale_factor = 2.0/largest_range;
        
        for(VertexIter vi = vertices.begin(); vi != vertices.end(); vi++){
            for(int i = 0; i < 3; i++){
                vi->coord[i] = (vi->coord[i] - center[i]) * scale_factor;
            }
        }
        
        
        return true;
    }
    
    template<> bool CMyMesh::ConstructMeshDataStructure(char *filename)
    {
        if( CMyMesh::Read(filename) == false ) return false;
        AddEdgeInfo();
    
        return true;
    }
    
    template<> void CMyMesh::Display(int mode)
    {
        glEnable(GL_LIGHTING);
        glEnable( GL_POLYGON_OFFSET_FILL );
        glPolygonOffset(1.0, 1.0);
        
        if(mode == 0){
            
            glEnable(GL_LIGHTING);
            
            for(FaceIter fi = faces.begin(); fi != faces.end(); fi++){
                
                if(fi->isActive){
                    glBegin(GL_TRIANGLES);
                    glNormal3dv(fi->halfedge[0].vertex->normal);
                    glVertex3dv(fi->halfedge[0].vertex->coord);
                    glNormal3dv(fi->halfedge[1].vertex->normal);
                    glVertex3dv(fi->halfedge[1].vertex->coord);
                    glNormal3dv(fi->halfedge[2].vertex->normal);
                    glVertex3dv(fi->halfedge[2].vertex->coord);
                    glEnd();
                }
            }
            
        }else if(mode == 1){
            
            for(FaceIter fi = faces.begin(); fi != faces.end(); fi++){
                
                if(fi->isActive){
                    glBegin(GL_TRIANGLES);
                    glNormal3dv(fi->halfedge[0].vertex->normal);
                    glVertex3dv(fi->halfedge[0].vertex->coord);
                    glNormal3dv(fi->halfedge[1].vertex->normal);
                    glVertex3dv(fi->halfedge[1].vertex->coord);
                    glNormal3dv(fi->halfedge[2].vertex->normal);
                    glVertex3dv(fi->halfedge[2].vertex->coord);
                    glEnd();
                }
            }
            
            glDisable(GL_LIGHTING);
            glLineWidth(1.0);
            glColor3d(0, 0, 0);
            
            for(FaceIter fi = faces.begin(); fi != faces.end(); fi++){
                if(fi->isActive){
                    glBegin(GL_LINE_LOOP);
                    glVertex3dv(fi->halfedge[0].vertex->coord);
                    glVertex3dv(fi->halfedge[1].vertex->coord);
                    glVertex3dv(fi->halfedge[2].vertex->coord);
                    glEnd();
                }
            }
        }
    }
}

#endif // !_MY_MESH_
