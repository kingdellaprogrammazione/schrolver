#include <vector>
#include <functional>
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <utility>
#include <gnuplot-iostream.h>
#include <fstream>
#include <string>

#define H_BAR 	1.054571817e-34 // J*s
#define M_STAR 9.1093837139e-31  // kg

// we use real numbers here -- can be demonstrated.

class Onedim_problem{
    double box_length; 
    int steps;      // insert even
    std::function<double(double)> potential;

    Eigen::VectorXd box_xvalues;
    Eigen::VectorXd psi;
    Eigen::VectorXd potential_values;
    Eigen::VectorXd eigenvalues; // dynamic vector

    Eigen::MatrixXd tridiagonal_matrix;  // dynamic matrix

    Gnuplot& gp;
    
    public:
    Onedim_problem(double sample_length,double sample_steps, std::function<double(double)> sample_potential, Gnuplot &sample_gp) 
    : box_length(sample_length), steps(sample_steps), potential(sample_potential), gp(sample_gp) {
        double step_length = box_length / steps;
        double center = box_length / 2;

        box_xvalues = Eigen::VectorXd::Zero(steps);

        for (int i = 0; i < steps; i++) {
            box_xvalues(i) = - center + i * step_length;
        } // first step is at the extremal

        tridiagonal_matrix = Eigen::MatrixXd::Zero(steps, steps);
        potential_values = Eigen::VectorXd::Zero(steps);
        psi = Eigen::VectorXd::Zero(steps);
        eigenvalues = Eigen::VectorXd::Zero(steps);

    }

    void discretize_potential(){
        for (int i =0 ; i < steps ; i++){
            double magnitude = potential(box_xvalues(i));
            potential_values(i) = magnitude;
        }
    }

    void build_matrix(){
        double step_length = box_length / steps;

        Eigen::VectorXd main_diagonal(steps);
        for (int i = 0; i < steps; i++) {
            main_diagonal(i) = H_BAR*H_BAR /( M_STAR * step_length * step_length) + 
            potential_values(i);
        }

        tridiagonal_matrix.diagonal() = main_diagonal;
        tridiagonal_matrix.diagonal(1) =  Eigen::VectorXd::Constant(steps-1, - H_BAR*H_BAR /( 2 * M_STAR * step_length * step_length) ); 
        tridiagonal_matrix.diagonal(-1) =  tridiagonal_matrix.diagonal(1);
        std::cout << tridiagonal_matrix << std::endl;
    }

    void compute_eigenvalues() {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(tridiagonal_matrix);
        eigenvalues = eigensolver.eigenvalues();
        std::cout << eigenvalues <<std::endl;
    }

};


double potential_well(double x, double half_well , double height) {
    if (abs(x)/(1./2) > half_well) {
        return height;
    } else {
        return 0.;
    }
}

double potential_well_proto(double x) {
    return potential_well(x, 0.1, 100 * 3.14 * 3.14 * H_BAR* H_BAR /(2* M_STAR *1* 1));
}

void print_with_gnu(){ // should be inside the onedim problem
    Gnuplot gp;
    
    // Define the range and step size for x values
    double x_min = -1.0e-10; // Start of the range (example)
    double x_max = 1.0e-10;  // End of the range (example)
    double step = 1.0e-12;   // Step size (example)
    
    // Generate x and y data points
    std::vector<std::pair<double, double>> points;
    for (double x = x_min; x <= x_max; x += step) {
        double y = potential_well_proto(x);
        points.emplace_back(x, y);
    }
    
    // Plot the function
    gp << "set title 'Potential Well Function'\n";
    gp << "set xlabel 'x'\n";
    gp << "set ylabel 'Potential'\n";
    gp << "plot '-' with lines title 'Potential Well'\n";
    gp.send1d(points);


}

void parse_input(std::string file_path, double &box_lenght, int &steps){
    std::ifstream file(file_path);  // Open the file for reading
    if (!file.is_open()) {
        std::cerr << "Failed to open the file.\n";
        exit(-1);
    }

    std::string line;
    std::getline(file, line);
        // Process each line here
        box_lenght = std::stod(line);

    std::getline(file, line);
        // Process each line here
    steps = std::stoi(line);

    file.close();  // Close the file

}


int main(void) {

    double box_lenght ;
    int steps ;
    // gp object should be defined here and passed to onedim_problem to be
    // able to print different onedim problem into the same window
    Gnuplot gp;

    parse_input("input.txt", box_lenght, steps);
    std::cout << box_lenght << steps << std::endl;
    
    Onedim_problem primo_problema = Onedim_problem(box_lenght, steps, potential_well_proto, gp);
    // primo_problema.discretize_potential();
    // primo_problema.build_matrix();
    //bprimo_problema.compute_eigenvalues();

    // to make gnuplot work, install boost library with boost devel, and 
    // compile linking boost lib: 
    // g++ 1_D.cpp -o main -lboost_iostreams -lboost_system -lboost_filesystem
    // do not use gnuplot-wx (buggy, perhaps https://sourceforge.net/p/gnuplot/bugs/2693/ related)
 


    return 0;
}

