#include "./Core/Solid.h"
#include "./Core/iterators.h"
#include "./Core/OBJFileReader.h"
#include "./Core/iterators.h"
#include "./Core/Point.h"
#include <math.h>
#include <fstream>
#include <vector>

using namespace MeshLib;

//  The following attributes have been added to the following classes

//  Vertex.h

    //  Point v_gradient;
    //  Point v_tgradient; // Gradient vector along the tangential space
    //  Point v_tgrad_f;
    //  Point v_search;
    //  Point v_normal;
    //  Point v_mp;

//  Face.h

    //  Point normal;
    //  double area;

//  Edge.h

    //  double kuv;


const double STEP_SIZE = 1e-2;
const double ENERGY_DIFF_THRESH = 1e-5;

Solid * mesh;

void setup(){

    // Setting up kuv

    for(SolidEdgeIterator eiter(mesh); !eiter.end(); ++eiter)
        (*eiter)->kuv = 1;
        
    for(SolidEdgeIterator eiter(mesh); !eiter.end(); ++eiter){

        // kuv = 0.5 (cot alpha + cot beta)
        // alpha and beta are the angles opposite to the edge

        Solid::tEdge edge = *eiter;

        Point p1 = mesh->edgeVertex1(edge)->point(); // edge point 1
        Point p2 = mesh->edgeVertex2(edge)->point(); // edge point 2
        Point p3 = edge->halfedge(0)->he_next()->target()->point(); // point 'towards' alpha
        Point p4 = edge->halfedge(0)->he_sym()->he_next()->target()->point(); // point 'towards' beta

        double alpha, beta;
        alpha = (p1 - p3)*(p2 - p3) / ((p1 - p3) ^ (p2 - p3)).norm() / 2; // 0.5 * cot alpha
        beta = (p1 - p4)*(p2 - p4) / ((p1 - p4) ^ (p2 - p4)).norm() / 2; // 0.5 * cot beta

        (*eiter)->kuv = alpha + beta;
    }

    // Setting up normals

	for (SolidFaceIterator fiter(mesh); !fiter.end(); ++fiter){

		Solid::tFace face = *fiter;
		Point faceVertex[3]; 

		int i = 0;
		for (FaceVertexIterator viter(face); !viter.end(); ++viter)
			faceVertex[i++] = (*viter)->point();

		Point normal = (faceVertex[1] - faceVertex[0]) ^ (faceVertex[2] - faceVertex[0]);
		face->normal = normal / normal.norm();
	}
}

void starMap(){

    // Calculating center point
    Point center(0, 0, 0);
    for(SolidVertexIterator viter(mesh); !viter.end(); ++viter){
        center += (*viter)->point();
    }
    center /= mesh->numVertices();

    // Normalizing all points
    for(SolidVertexIterator viter(mesh); !viter.end(); ++viter){
		(*viter)->point() -= center;
        (*viter)->point() /= (*viter)->point().norm();
        (*viter)->v_normal = (*viter)->point();
	}
}

double getEnergy(bool HARMONIC = true){

    double energy{0};
    for(SolidEdgeIterator eiter(mesh); !eiter.end(); ++eiter){

        Solid::tVertex v1 = mesh->edgeVertex1(*eiter);
        Solid::tVertex v2 = mesh->edgeVertex2(*eiter);
        Point uv = v1->point() - v2->point();

        if(HARMONIC){
            // std::cout << "getEnergy HARMONIC, energy += " << (*eiter)->kuv * uv.norm2() << std::endl;
            energy += (*eiter)->kuv * uv.norm2();
        }
        else{
            // std::cout << "getEnergy TUETTE, energy += " << (*eiter)->kuv * uv.norm2() << std::endl;
            energy += uv.norm2();
        }
    }
    return energy;
}

void setupGradient(bool HARMONIC = true){

    for(SolidVertexIterator viter(mesh); !viter.end(); ++viter){

        Solid::tVertex i = *viter;
        Point gradient(0, 0, 0);
        for(VertexVertexIterator vviter(i); !vviter.end(); ++vviter){
            Vertex *j = *vviter;
            if(HARMONIC)
                gradient += (i->point() - j->point()) * mesh->vertexEdge(i, j)->kuv;
            else
                gradient += (i->point() - j->point());
        }
        (*viter)->v_gradient = gradient;
        (*viter)->v_tgradient = gradient - (*viter)->v_normal * (gradient * ((*viter)->v_normal));
    }
}

void updateMesh(){

    for(SolidVertexIterator viter(mesh); !viter.end(); ++viter){
        Solid::tVertex v = *viter;
        v->point() -= (*viter)->v_tgradient * STEP_SIZE;
        v->point() /= v->point().norm();
    }
}

void updateCOM(){

    for(SolidFaceIterator fiter(mesh); !fiter.end(); ++fiter){

        Solid::tFace face = *fiter;
        Point faceVertex[3];

        int i{0};
        for(FaceVertexIterator viter(face); !viter.end(); ++viter, ++i)
            faceVertex[i] = (*viter)->point();

        Point n = (faceVertex[1] - faceVertex[0]) ^ (faceVertex[2] - faceVertex[0]);
        face->area = n.norm() / 2.0;
    }
    for (SolidVertexIterator viter(mesh); !viter.end(); ++viter){

		Vertex *vertex = *viter;
		double area{0};

		for (VertexFaceIterator fiter(vertex); !fiter.end(); ++fiter)
			area += (*fiter)->area;

		(*viter)->v_area = area / 3.0;
	}

	Point center(0, 0, 0);
	double mass = 0;

	for (SolidVertexIterator viter(mesh); !viter.end(); ++viter) {
		Vertex *vertex = *viter;
		center += vertex->point() * (*viter)->v_area;
		mass += (*viter)->v_area;
	}
	center /= mass;

	for (SolidVertexIterator viter(mesh); !viter.end(); ++viter) {
		Vertex *vertex = *viter;
		vertex->point() -= center;
		vertex->point() /= vertex->point().norm();
	}
}

void conjugateGradientMethod(){

    for(SolidVertexIterator viter(mesh); !viter.end(); ++viter){
        (*viter)->v_beta = ((*viter)->v_tgradient * (*viter)->v_tgradient) / ((*viter)->v_tgrad_f * (*viter)->v_tgrad_f);
        (*viter)->v_search = (*viter)->v_search * (*viter)->v_beta - (*viter)->v_tgradient;
    }
    double alpha = 1e-6;
    double init_energy = getEnergy(true);

    for(SolidVertexIterator viter(mesh); !viter.end(); ++viter)
        (*viter)->v_mp = (*viter)->point();

    while(alpha < 1e-2){

        for(SolidVertexIterator viter(mesh); !viter.end(); ++viter){
            Solid::tVertex v = *viter;
            v->point() += (*viter)->v_search * alpha;
            v->point() /= v->point().norm();
        }
        double cur_energy = getEnergy(true);
        if(cur_energy < init_energy){
            for(SolidVertexIterator viter(mesh); !viter.end(); ++viter)
                (*viter)->v_mp = (*viter)->point();
        }

        alpha *= 2;
    }

    for(SolidVertexIterator viter(mesh); !viter.end(); ++viter)
        (*viter)->point() = (*viter)->v_mp;
}

void tuetteMap(){

    double te = getEnergy(false); // tuette energy
    double oe = 1000;

    for(int iterations{0}; oe-te > ENERGY_DIFF_THRESH; iterations++){
        setupGradient(false);
		updateMesh();
		oe = te;
		te = getEnergy(false);
		if (iterations % 100 == 0 && iterations < 1000 || iterations % 1000 == 0){
            std::cout << "Iteration number: " << iterations << std::endl;
            std::cout << "Energy: " << te << std::endl;
            std::cout << "Energy difference: " << oe - te << std::endl;
            std::cout << std::endl;
        }
    }
}

void harmonicMap(){

    double te = getEnergy(true); // tuette energy
    double oe = 1000;

    setupGradient(true);
    updateMesh();

    for(SolidVertexIterator viter(mesh); !viter.end(); ++viter){
        (*viter)->v_tgrad_f = (*viter)->v_tgradient;
        (*viter)->v_search = (*viter)->v_tgradient;
    }
    double he = getEnergy(true); // harmonic energy

    for(int iterations{0}; oe - he > ENERGY_DIFF_THRESH; iterations++){
        setupGradient(true);
        conjugateGradientMethod();
        oe = he;
        he = getEnergy(true);

        if(iterations < 10 or iterations % 100 == 0){
            std::cout << "Iteration number: " << iterations << std::endl;
            std::cout << "Energy: " << he << std::endl;
            std::cout << "Energy difference: " << oe - he << std::endl;
            std::cout << std::endl;
        }
    }
}

int main(int argc, char * argv[]){

    // assert(0);

    Solid inp_mesh;
	OBJFileReader of;
	std::ifstream in(argv[1]);
	of.readToSolid(&inp_mesh, in);

    mesh = &inp_mesh;

    setup();
    starMap();

    std::string out;

    out = std::string(argv[2]) + "Star.obj";
    of.writeToObj(mesh, out);

    std::cout << "\nGenerating Tuette Map" << std::endl;

    tuetteMap();
    out = std::string(argv[2]) + "Tuette.obj";
	of.writeToObj(mesh, out);

    std::cout << "\n\nGenerating Harmonic Map" << std::endl;

    harmonicMap();
    out = std::string(argv[2]) + "Harmonic.obj";
	of.writeToObj(mesh, out);

    return 0;
}