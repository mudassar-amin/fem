#include "F:/Course/fem/eigen/Eigen/Eigen"
#include "F:/Course/fem/eigen/Eigen/Dense"
#include <ttl/halfedge/HeTriang.h>
#include <ttl/halfedge/HeDart.h>
#include <ttl/halfedge/HeTraits.h>
#define _USE_MATH_DEFINES
#include <fstream>
#include <iostream>
#include <iomanip> // For std::setw

#include <algorithm>
#include <math.h>


using namespace std;
using namespace hed;
using namespace Eigen;



class FEMobject {
private:
    vector<Node*>* nodes = new vector<Node*>;
    Triangulation triang;
    Eigen::SparseMatrix<double> A;
    Eigen::SparseMatrix<double> M;
    Eigen::SparseMatrix<double> R;
    Eigen::VectorXd b;
    Eigen::VectorXd r;

public:
    std::string problemtype;

    void CircleMesh(int n, int m, double r){
           Node* p = new Node(0, 0);
           nodes->push_back(p);

           for (int k = 0;  k < n; k++){
               for (int i = 0; i <(m *(k + 1)); i++ ){
                   auto vec  = Eigen::Vector2d((r *( k + 1))/n , 0);
                   auto a = (i * 2 * M_PI)/(m * (k + 1));

                   Eigen::Matrix2d v;
                   v << cos(a), -sin(a),
                        sin(a),  cos(a);

                   auto t = v * vec;

                   p = new Node(t.x(), t.y());
                   nodes->push_back(p);

               }
           }
            triang.createDelaunay(nodes->begin(), nodes->end());
       }


    void SquareMesh(int n, double d, Eigen::Vector2d Op) {
        for (int j = 0; j <= n; ++j) {
            for (int i = 0; i <= n; ++i) {
                double x = Op(0) + ((i * d) / n); // x-coordinate
                double y = Op(1) + ((j * d) / n); // y-coordinate
                hed::Node* p = new hed::Node(x, y);
                // Add the new node to the nodes vector
                nodes->push_back(p);
            }
        }
        triang.createDelaunay(nodes->begin(), nodes->end());
    }





    double kappa(double x, double y){
        if(this->problemtype == "laplace") {
            return 1e6;
        } else if(this->problemtype == "poisson") {
            return 1e6;
        }else if(this->problemtype == "helmholtz") {
            if (x > 0.0)
                return 0.0;
            else
                return 1e6;
        }
        return 0; //default value
    }

    double gN(double x, double y) {
        if(this->problemtype == "laplace") {
            return 0.0;
        } else if(this->problemtype == "poisson") {
            return 0.0;
        }else if(this->problemtype == "helmholtz") {
            return 0.0;
        }
        return 0; //default value
        }

    double gD(double x, double y) {
        if (this->problemtype == "laplace") {
            double phi = atan2(y, x);
            return cos(4 * phi);
        } else if (this->problemtype == "poisson") {
            return y * y / 2.0; // y squared divided by 2 for Poisson's problem
        } else if (this->problemtype == "helmholtz") {
            return 0.25; // Constant boundary value for Helmholtz problem
        }
        return 0;
    }

    double f(double x, double y) {

            return 1.0;
        }

    double triarea(Node* N1, Node* N2, Node* N3) {
        Vector3d a = {N2->x() - N1->x(), N2->y() - N1->y(), 0};
        Vector3d b = {N3->x() - N1->x(), N3->y() - N1->y(), 0};
        double area = 0.5*(b.cross(a)).norm();
        return area;
        }

    Vector2<Vector3d> gradients(Vector3d x, Vector3d y, double area) {
        Vector3d b((y(1)-y(2))/(2*area), (y(2)-y(0))/(2*area), (y(0)-y(1))/(2*area));
        Vector3d c((x(2)-x(1))/(2*area), (x(0)-x(2))/(2*area), (x(1)-x(0))/(2*area));
        Vector2<Vector3d> gradphi = {b, c};
        return gradphi;

        }


    SparseMatrix<double> stiffMat(list<Edge*> trilist, int np) {
        SparseMatrix<double> A(np, np);
        list<Edge*>::iterator K;
        for (K = trilist.begin(); K != trilist.end(); K++) {
            Edge* edg = *K; // extract the edge of the triangle
            // extract nodes...
            Node* N1 = edg->getSourceNode();
            Node* N2 = edg->getTargetNode();
            Node* N3 = edg->getNextEdgeInFace()->getTargetNode();
            // ...and their global indices
            Vector3<Index> loc2glb = {N1->id(), N2->id(), N3->id()};
            // extract coordinates of the nodes
            Vector3d x = {N1->x(), N2->x(), N3->x()};
            Vector3d y = {N1->y(), N2->y(), N3->y()};
            double area = triarea(N1, N2, N3); // compute the triangle area
            Vector2<Vector3d> gradphi = gradients(x, y, area); // compute the gradients
            // compute AK
            Vector3d b = gradphi(0);
            Vector3d c = gradphi(1);
            Matrix3d AK = (b * b.transpose() + c * c.transpose()) * area;
            // insert AK to the global stiffness matrix A
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    A.coeffRef(loc2glb(i), loc2glb(j)) += AK(i, j);
                }
            }
        }
        return A; // Move this line outside the loop
    }


    SparseMatrix<double> massMat(list<Edge*> trilist, int np) {
        SparseMatrix<double> M(np, np); // initialization of the mass matrix M
        list<Edge*>::iterator K;
        for (K = trilist.begin(); K != trilist.end(); K++) { // loop over triangles
            Edge* edg = *K; // extract the edge of the triangle
            // extract nodes...
            Node* N1 = edg->getSourceNode();
            Node* N2 = edg->getTargetNode();
            Node* N3 = edg->getNextEdgeInFace()->getTargetNode();
            // ...and their global indices
            Vector3<Index> loc2glb = {N1->id(), N2->id(), N3->id()};
            // extract coordinates of the nodes
            Vector3d x = {N1->x(), N2->x(), N3->x()};
            Vector3d y = {N1->y(), N2->y(), N3->y()};
            double area = triarea(N1, N2, N3); // compute the triangle area
            // compute MK
            Matrix3d mat;
            mat << 2, 1, 1,
                   1, 2, 1,
                   1, 1, 2;
            Matrix3d MK = (1.0 / 12) * mat * area;
            // insert MK to the global mass matrix M
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    M.coeffRef(loc2glb(i), loc2glb(j)) += MK(i, j);
                }
            }
        }
        return M;
    }


    VectorXd loadVect(list<Edge*> trilist, int np) {
        VectorXd b = VectorXd::Zero(np); // initialization of the load vector b
        list<Edge*>::iterator K;
        for (K = trilist.begin(); K != trilist.end(); K++) { // loop over triangles
            Edge* edg = *K; // extract the edge of the triangle
            // extract nodes...
            Node* N1 = edg->getSourceNode();
            Node* N2 = edg->getTargetNode();
            Node* N3 = edg->getNextEdgeInFace()->getTargetNode();
            // ...and their global indices
            Vector3<Index> loc2glb = {N1->id(), N2->id(), N3->id()};
            // extract coordinates of the nodes
            Vector3d x = {N1->x(), N2->x(), N3->x()};
            Vector3d y = {N1->y(), N2->y(), N3->y()};
            double area = triarea(N1, N2, N3); // compute the triangle area
            // compute bK
            Vector3d F = {f(x(0), y(0)), f(x(1), y(1)), f(x(2), y(2))};
            Vector3d bK = F / 3 * area;
            // add bK to the global load vector b
            for (int j = 0; j < 3; ++j) {
                b(loc2glb(j)) += bK(j);
            }
        }
        return b;
    }


    SparseMatrix<double> RobinMat(list<Dart> boundary, int np) {
        SparseMatrix<double> R(np, np); // initialization of the boundary matrix R
        list<Dart>::iterator E;
        for (E = boundary.begin(); E != boundary.end(); ++E) { // loop over boundary edges
            // extract nodes...
            Edge* edg = E->getEdge();
            Node* N1 = edg->getSourceNode();
            Node* N2 = edg->getTargetNode();
            // ...and their global indices
            Vector2<Index> loc2glb = {N1->id(), N2->id()};
            // extract coordinates of the nodes
            Vector2d x = {N1->x(), N2->x()};
            Vector2d y = {N1->y(), N2->y()};
            // compute the length of the edge
            double len = sqrt((x(0) - x(1)) * (x(0) - x(1)) + (y(0) - y(1)) * (y(0) - y(1)));
            // find edge centroid
            double xc = (x(0) + x(1)) / 2;
            double yc = (y(0) + y(1)) / 2;
            double k = kappa(xc, yc); // compute the value of κ at centroid
            // compute RE
            Matrix2d mat;
            mat << 2, 1,
                1, 2;
            Matrix2d RE = k / 6 * mat * len;
            // add RE to the global boundary matrix R
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    R.coeffRef(loc2glb(i), loc2glb(j)) += RE(i, j);
                }
            }
        }
        return R;
    }


    VectorXd RobinVect(list<Dart> boundary, int np) {
        VectorXd r = VectorXd::Zero(np); // initialization of the boundary vector r
        list<Dart>::iterator E;
        for (E = boundary.begin(); E != boundary.end(); ++E) { // loop over boundary edges

            Edge* edg = E->getEdge();
            Node* N1 = edg->getSourceNode();
            Node* N2 = edg->getTargetNode();

            Vector2<Index> loc2glb = {N1->id(), N2->id()};

            Vector2d x = {N1->x(), N2->x()};
            Vector2d y = {N1->y(), N2->y()};

            double len = sqrt((x(0) - x(1)) * (x(0) - x(1)) + (y(0) - y(1)) * (y(0) - y(1)));

            double xc = (x(0) + x(1)) / 2;
            double yc = (y(0) + y(1)) / 2;
            double tmp = kappa(xc, yc) * gD(xc, yc) + gN(xc, yc);
            Vector2d vec = {1, 1};
            Vector2d rE = tmp * vec * len / 2;

            for (int j = 0; j < 2; ++j) {
                r(loc2glb(j)) += rE(j);
            }
        }
        return r;
    }



    void solve() {
        int np = nodes->size();
        auto trilist = triang.getLeadingEdges();
        auto nodelist = triang.getNodes();
        Edge* edge = triang.getBoundaryEdge();
        Dart b_dart(edge); //creates a boundary dart
        list<Dart> boundary;
        ttl::getBoundary(b_dart, boundary);

        // Initialize the solver
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        Eigen::VectorXd zeta; // Solution vector

        // Assemble the stiffness matrix A and Robin matrix R
        A = stiffMat(trilist, np);
        R = RobinMat(boundary, np);
        b = loadVect(trilist, np);
        r = RobinVect(boundary, np);

        if (this->problemtype == "laplace") {
           Eigen::VectorXd b = Eigen::VectorXd::Zero(np);
            zeta = solver.compute(A + R).solve(r);

        } else if (this->problemtype == "poisson") {
            Eigen::VectorXd b = loadVect(trilist, np);
            zeta = solver.compute(A + R).solve(r+b);

        } else if (this->problemtype == "helmholtz") {

            double lambda = 81;
            M = massMat(trilist, np);
            Eigen::VectorXd r = RobinVect(boundary, np);
            zeta = solver.compute(A + R - (lambda * M)).solve(r);
        } else {
            std::cerr << "Unknown problem type: " << this->problemtype << std::endl;
            return;
        }

        //cout<<zeta;

        // Update node values
        list<Node*>::iterator L;
        for (L = nodelist->begin(); L != nodelist->end(); L++) {
            Node* node = *L;
            node->init(node->x(),node->y(), zeta(node->id())); //set the node values
        }

    }


    void visualization(const std::string& filename) {
            std::ofstream objfile(filename);
            if (!objfile.is_open()) {
                std::cerr << "Failed to open " << filename << " for writing.\n";
                return;
            }

            // Get a list of nodes and edges from the triangulation
            auto nodelist = triang.getNodes();
            auto trilist = triang.getLeadingEdges();
            int np = nodes->size();

            // Save vertex positions in the OBJ file
            int normind = 1;
            for (const auto& node : *nodelist) {
                objfile << "v " << node->x() << " " << node->z() <<  " " << node->y()<< "\n"; // Z-coordinate is set to 0 for 2D
            }  //-uexact(node->x(), node->y())


            objfile << "vn 0.0 0.0 1.0\n";

            // Save faces
            for (const auto& edge : trilist) {

                auto N1 = edge->getSourceNode();
                auto N2 = edge->getTargetNode();
                auto N3 = edge->getNextEdgeInFace()->getTargetNode();

                objfile << "f " << N1->id()-np-3 << "//" << normind << " " <<
                N2->id()-np-3 << "//" << normind << " " <<
                N3->id()-np-3 << "//" << normind << "\n"; // f v1//vn1 v2//vn2 v3//vn3
                normind++;
            }

            objfile.close();
            if (!objfile.good()) {
                std::cerr << "An error occurred while writing to " << filename << "\n";
            }
        }



    double uexact(double x, double y){

        if(this->problemtype == "laplace") {
            double rho = sqrt((x * x) + (y * y));
            double phi = atan2(y, x);
            return pow(rho, 4) * cos(4 * phi);

        } else if(this->problemtype == "poisson") {
            return (1 - (x * x)) / 2.0;

        }else if(this->problemtype == "helmholtz") {
            const double lambda = 81.0;
            const double sqrt_lambda = sqrt(lambda);
            return (cos(sqrt_lambda * x) + tan(sqrt_lambda) * sin(sqrt_lambda * x))*0.25;

        }
        return 0; //default value
    }

    double getError(){
        double error = 0;
        list<Edge*>::iterator K;
        //auto nodelist = triang.getNodes();
        auto trilist = triang.getLeadingEdges();
        // Assuming you have a collection of edges to iterate over
        for (K = trilist.begin(); K != trilist.end(); K++) {
        Edge* edge = *K;
        Node* nodea = edge->getSourceNode();
        Node* nodeb = edge->getTargetNode();
        Node* nodec = edge->getNextEdgeInFace()->getTargetNode();
        Eigen::Vector3d x = {nodea->x(), nodeb->x(), nodec->x()};
        Eigen::Vector3d y = {nodea->y(), nodeb->y(), nodec->y()};
        double area = triarea(nodea, nodeb, nodec);
        double xc = (x(0) + x(1) + x(2)) / 3;
        double yc = (y(0) + y(1) + y(2)) / 3;
        double ubar = uexact(xc, yc);
        double uhbar = (nodea->z() + nodeb->z() + nodec->z()) / 3;
        double eK = (ubar - uhbar) * (ubar - uhbar) * area; // Corrected exponentiation
        error += eK;
        }

        return sqrt(error);
    }

    int getDoFs(){

       return this->triang.getNodes()->size(); // This is correct
    }

};


int main() {

    int no_of_nodes =25;
    Eigen::Vector2d Op(0, -0.5);
    FEMobject model;

    //model.problemtype = "laplace";
    //model.problemtype = "poisson";
    model.problemtype = "helmholtz";


    model.SquareMesh(no_of_nodes, 1, Op );
    //model.CircleMesh(5, no_of_nodes, 1);
    model.solve();


    //model.visualization("squarelaplace.obj");
    //model.visualization("circlelaplace.obj");


    //model.visualization("squarepoisson.obj");
    //model.visualization("circlepoisson.obj");


    model.visualization("squarehelmholtz.obj");
    //model.visualization("circlehelmholtz.obj");






    int numDoFs = model.getDoFs();
        double error = model.getError();
        std::string problemType = model.problemtype; // Ensure this is set correctly in FEMobject

            // Print the header of the table
            std::cout << std::left << std::setw(15) << "Problem Type"
                      << std::setw(15) << "DoFs"
                      << std::setw(15) << "Error" << std::endl;

            // Print a line for the table
            std::cout << std::setfill('-') << std::setw(45) << "" << std::endl;
            std::cout << std::setfill(' '); // Reset the fill character

            // Print the data in a tabular format
            std::cout << std::left << std::setw(15) << model.problemtype
                      << std::setw(15) << numDoFs
                      << std::setw(15) << std::setprecision(5) << error << std::endl;



  return 0;
}






