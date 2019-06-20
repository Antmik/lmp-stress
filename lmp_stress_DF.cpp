//Examples:
//lmp_stress filename dump.atom NVT n_steps 10 n_bins 100 kernel gaussian 0.1 density momentum viscosity
//lmp_stress filename dump.atom NVT n_steps 10 n_bins 100 kernel gaussian 0.1 all

#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include <math.h>
#include <time.h>

double cubic_bond_function(const double pos, const double R_cut ) {
    if (pos >= 0.0 && pos<R_cut) {
        return ( pos - pos*pos*pos/R_cut/R_cut + pos*pos*pos*pos/2/R_cut/R_cut/R_cut );
    }
    else if (pos >= 0.0 && pos> R_cut) {
        return R_cut/2;
    }
    else if (pos < 0.0 && pos> -R_cut) {
        return -( abs(pos) - abs(pos*pos*pos)/R_cut/R_cut + abs(pos*pos*pos*pos)/2/R_cut/R_cut/R_cut);
    }
    else if (pos < 0.0 && pos< -R_cut) {
        return -R_cut/2;
    }
    else{
        return 0.0 ;
    }
}


double kernel(arma::vec& dist, arma::vec& Prefactor, const double R_cut, const double flag_kernel) {
	double result;
	if (flag_kernel == 0) {
		result = arma::sum(0.5/(2.0* R_cut)*Prefactor % ( arma::sign(dist + R_cut)- arma::sign(dist - R_cut)));
			//exp(-(dist_sq) / (2 * R_cut*R_cut)) / sqrt(2 * arma::datum::pi *R_cut*R_cut));
	}
	else if (flag_kernel==1) {
//        arma::vec dist_sq = dist % dist;
		result = arma::sum(Prefactor % arma::exp(-(dist % dist) / (2 * R_cut*R_cut)) / sqrt(2 * arma::datum::pi *R_cut*R_cut));
	}
	else if (flag_kernel==2) {
        arma::vec dist_tmp = arma::abs(dist);
        
        arma::vec result_tmp= 1/R_cut- (3/(R_cut*R_cut*R_cut))*( dist_tmp % dist_tmp ) + (2/(R_cut*R_cut*R_cut*R_cut))*( (dist_tmp % dist_tmp) % dist_tmp);
        
        for (int i = 0; i < dist_tmp.n_elem; i++) {
            if (dist_tmp(i) > R_cut){
            result_tmp(i) = 0.0;
            }
        }
        result = arma::sum( Prefactor % result_tmp) ;
	}

	return result;
}

double kernel_derivative(arma::vec& dist, arma::vec& Prefactor, const double R_cut, const double flag_kernel) {
    double result;
    if (flag_kernel == 0) {
        std::cout << "No kernel derivative for constant type";
    }
    else if (flag_kernel==1) {
        arma::vec dist_sq = dist % dist;
        result = arma::sum(Prefactor % (- dist / (R_cut*R_cut) ) % arma::exp(-(dist_sq) / (2 * R_cut*R_cut)) / sqrt(2 * arma::datum::pi *R_cut*R_cut));
    }
    else if (flag_kernel==2) {
        arma::vec dist_tmp = dist;
        arma::vec result_tmp=dist;
        
        for (int i = 0; i < dist_tmp.n_elem; i++) {
            if (dist_tmp(i) > R_cut){
                result_tmp(i) = 0.0;
            }
            if (dist_tmp(i) >= 0.0 && dist_tmp(i)<R_cut) {
                result_tmp(i) = - 3*2*dist_tmp(i)/R_cut/R_cut/R_cut + 2*3*dist_tmp(i)*dist_tmp(i)/R_cut/R_cut/R_cut/R_cut;
            }
            else if (dist_tmp(i) >= 0.0 && dist_tmp(i)> R_cut) {
                result_tmp(i) = 0.0;
            }
            else if (dist_tmp(i) < 0.0 && dist_tmp(i)> -R_cut) {
               result_tmp(i) = - ( 3*2*dist_tmp(i)/R_cut/R_cut/R_cut + 2*3*dist_tmp(i)*dist_tmp(i)/R_cut/R_cut/R_cut/R_cut);
            }
            else if (dist_tmp(i) < 0.0 && dist_tmp(i)< -R_cut) {
                result_tmp(i) = 0.0;
            }
            else{
                result_tmp(i) = 0.0;
            }
        }
        
        result = arma::sum( Prefactor % result_tmp) ;
    }
    
    return result;
}

int compute_density(const int n_steps, const int n_bin, const double Lx, const double Lz, const double R_cut, arma::vec& Mass, arma::vec& density, arma::mat& y_vec, arma::vec& bin_vec, const double flag_kernel) {
	for (int i_bin = 0; i_bin < n_bin; i_bin++) {
		for (int i_step = 0; i_step < n_steps; i_step++) {
			arma::vec dist = bin_vec(i_bin) - y_vec.col(i_step);
			density(i_bin) += (1 / Lx / Lz * kernel ( dist , Mass, R_cut, flag_kernel)) / n_steps;
	}
	}
	return 0;
}

int compute_density_derivative(const int n_steps, const int n_bin, const double Lx, const double Lz, const double R_cut, arma::vec& Mass, arma::vec& density_derivative, arma::mat& y, arma::vec& bin_vec, const double flag_kernel) {
    for (int i_bin = 0; i_bin < n_bin; i_bin++) {
        for (int i_step = 0; i_step < n_steps; i_step++) {
            arma::vec dist = bin_vec(i_bin) - y.col(i_step);
            density_derivative(i_bin) += (1 / Lx / Lz * kernel_derivative ( dist , Mass, R_cut, flag_kernel)) / n_steps;
        }
    }
    return 0;
}

int compute_momentum(const int n_steps, const int n_bin, const double Lx, const double Lz, const double R_cut, arma::mat& Velx, arma::vec& momentum, arma::mat& y, arma::vec& bin_vec, const double flag_kernel) {
	for (int i_bin = 0; i_bin < n_bin; i_bin ++) {
		for (int i_step = 0; i_step < n_steps; i_step++) {
			arma::vec dist = bin_vec(i_bin) - y.col(i_step);
			arma::vec vel = Velx.col(i_step);
			momentum(i_bin) += (1 / Lx / Lz * kernel(dist, vel, R_cut, flag_kernel)) / n_steps;
		}
	}
	return 0;
}

int compute_momentum_y(const int n_steps, const int n_bin, const double Lx, const double Lz, const double R_cut, arma::mat& Vely, arma::vec& momentum_y, arma::mat& y, arma::vec& bin_vec, const double flag_kernel) {
    for (int i_bin = 0; i_bin < n_bin; i_bin ++) {
        for (int i_step = 0; i_step < n_steps; i_step++) {
            arma::vec dist = bin_vec(i_bin) - y.col(i_step);
            arma::vec vel = Vely.col(i_step);
            momentum_y(i_bin) += (1 / Lx / Lz * kernel(dist, vel, R_cut, flag_kernel)) / n_steps;
        }
    }
    return 0;
}

int compute_momentum_derivative(const int n_steps, const int n_bin, const double Lx, const double Lz, const double R_cut, arma::mat& Velx, arma::vec& momentum_derivative, arma::mat& y, arma::vec& bin_vec, const double flag_kernel) {
    for (int i_bin = 0; i_bin < n_bin; i_bin ++) {
        for (int i_step = 0; i_step < n_steps; i_step++) {
            arma::vec dist = bin_vec(i_bin) - y.col(i_step);
            arma::vec vel = Velx.col(i_step);
            momentum_derivative(i_bin) += (1 / Lx / Lz * kernel_derivative(dist, vel, R_cut, flag_kernel)) / n_steps;
        }
    }
    return 0;
}

int compute_momentum_y_derivative(const int n_steps, const int n_bin, const double Lx, const double Lz, const double R_cut, arma::mat& Vely, arma::vec& momentum_y_derivative, arma::mat& y, arma::vec& bin_vec, const double flag_kernel) {
    for (int i_bin = 0; i_bin < n_bin; i_bin ++) {
        for (int i_step = 0; i_step < n_steps; i_step++) {
            arma::vec dist = bin_vec(i_bin) - y.col(i_step);
            arma::vec vel = Vely.col(i_step);
            momentum_y_derivative(i_bin) += (1 / Lx / Lz * kernel_derivative(dist, vel, R_cut, flag_kernel)) / n_steps;
        }
    }
    return 0;
}

int compute_temperature(const int n_steps, const int n_bin, const double Lx, const double Lz, const double R_cut, arma::mat& velx, arma::mat& vely, arma::mat& velz, arma::mat& ave_velx,arma::mat& ave_vely, arma::vec& temp, arma::mat& y, arma::vec& bin_vec, const double flag_kernel) {
	for (int i_bin = 0; i_bin < n_bin; i_bin++) {
		for (int i_step = 0; i_step < n_steps; i_step++) {
			arma::vec dist = bin_vec(i_bin) - y.col(i_step);
			arma::vec vel = (velx.col(i_step) - ave_velx(i_bin)) % (velx.col(i_step) - ave_velx(i_bin)) + (vely.col(i_step) - ave_vely(i_bin)) % (vely.col(i_step) - ave_vely(i_bin)) + (velz.col(i_step) % velz.col(i_step));
			temp(i_bin) += (1 / 3.0 / Lx / Lz * kernel(dist, vel, R_cut, flag_kernel)) / n_steps;
		}
	}
	return 0;
}

int compute_stressK(const int n_steps, const int n_bin, const double Lx, const double Lz, const double R_cut, arma::mat& vel1, arma::mat& vel2, arma::mat& ave_vel1, arma::mat& ave_vel2, arma::vec& stressK, arma::mat& y, arma::vec& bin_vec, const double flag_kernel) {
	for (int i_bin = 0; i_bin < n_bin; i_bin++) {
		for (int i_step = 0; i_step < n_steps; i_step++) {
			arma::vec dist = bin_vec(i_bin) - y.col(i_step);
			arma::vec vel = (vel1.col(i_step) - ave_vel1(i_bin)) % (vel2.col(i_step) - ave_vel2(i_bin));
			stressK(i_bin) -= (1 / Lx / Lz * kernel(dist, vel, R_cut, flag_kernel)) / n_steps;
		}
	}
	return 0;
}

int read_file(const std::string& namefile, const int n_steps, const int n_particles, arma::mat& type , arma::mat& x, arma::mat& y, arma::mat& z, arma::mat& vx, arma::mat& vy, arma::mat& vz) {
	std::string line;
	//ifstream myfile("example.txt");
	std::ifstream myfile(namefile);
	if (myfile.is_open())
	{
		int n_line = 0;
		int n_time = 0;
		while (std::getline(myfile, line) && n_time<n_steps)
		{
			if (n_line < 9) {
				n_line += 1;
			}
			else if (n_line == n_particles + 8) {
				n_time += 1;
				n_line = 0;
			}
			else {

				//split line
				std::vector<std::string> result;
				std::istringstream iss(line);
				for (std::string line; iss >> line;) {
					result.push_back(line);
				}

				type(n_line - 8, n_time) = stod(result[1])-1;
				x(n_line - 8, n_time) =	stod(result[2]);
				y(n_line - 8, n_time) = stod(result[3]);
				z(n_line - 8, n_time) = stod(result[4]);
				vx(n_line - 8, n_time) = stod(result[5]);
				vy(n_line - 8, n_time) = stod(result[6]);
				vz(n_line - 8, n_time) = stod(result[7]);
				n_line += 1;
				//std::cout << line << '\n';
			
			}
		}
		myfile.close();
	}
	else std::cout << "Unable to open file";
	return 0;
}

int read_parameters(const std::string& namefile, int &n_types, int &n_particles, double &Lx, double &Ly, double &Lz) {
	std::string line;
	std::ifstream myfile(namefile);
	if (myfile.is_open())
	{
		int n_line = 0;
		n_particles = 0;
		n_types = 1;
		while (std::getline(myfile, line) && n_line< n_particles+8)
		{
			if (n_line == 3) {
				n_particles = stoi(line);
				n_line += 1;
			}
			else if (n_line == 5) {
				std::vector<std::string> result;
				std::istringstream iss(line);
				for (std::string line; iss >> line;) {
					result.push_back(line);
				}
				Lx = stod(result[1]) - stod(result[0]);
				n_line += 1;
			}
			else if (n_line == 6) {
				std::vector<std::string> result;
				std::istringstream iss(line);
				for (std::string line; iss >> line;) {
					result.push_back(line);
				}
				Ly = stod(result[1]) - stod(result[0]);
				n_line += 1;
			}
			else if (n_line == 7) {
				std::vector<std::string> result;
				std::istringstream iss(line);
				for (std::string line; iss >> line;) {
					result.push_back(line);
				}
				Lz = stod(result[1]) - stod(result[0]);
				n_line += 1;
			}
			else if (n_line > 9) {
				std::vector<std::string> result;
				std::istringstream iss(line);
				for (std::string line; iss >> line;) {
					result.push_back(line);
				}
				if (stod(result[1]) > n_types) {
					n_types = stoi(result[1]);
				}
				n_line += 1;
			}
			else {
				n_line += 1;
			}
		}
		myfile.close();
	}
	else std::cout << "Unable to open file";
	return 0;
	}


int read_parameter_file(const std::string& namefile,  arma::mat& epsilon , double &Force_cut) {
	std::string line;
	std::ifstream myfile(namefile);
	if (myfile.is_open())
	{
		while (std::getline(myfile, line) )
		{
				//split line
				std::vector<std::string> result;
				std::istringstream iss(line);
				for (std::string line; iss >> line;) {
					result.push_back(line);
				}
				
				if (&result[0] != NULL){
					std::string str1= result[0];
					//std::cout << "str1" << str1 << "\n";
						if (str1.compare("pair_coeff")==0){
						std::string stri= result[1];
						std::string strj= result[2];
						if (stri.compare("*")==0 || strj.compare("*")==0){
							Force_cut=	stod(result[5]);
							std::cout << "All epsilon set to 1" << '\n';
						}
						else
						{
							int i=stoi(result[1])-1;
							int j=stoi(result[2])-1;
							double e_value=stod(result[3]);
							epsilon(i,j)=e_value;
							epsilon(j,i)=e_value;
							Force_cut=	stod(result[5]);
						}
					}
				}

		}
		myfile.close();
		std::cout << "Epsilon " << epsilon << "\n";
		std::cout << "Force_cut " << Force_cut << "\n";
	}
	else{
		std::cout << "Unable to open file";
	}
	return 0;
}

int main(int argc, const char **argv) {
	clock_t tStart = clock();
	
	std::string filename;
	std::string parameter_filename;
	int	n_steps;
	int	n_bin;

	int n_types=0;
	int n_particles=0;
	double	Lx=0;
	double	Ly=0;
	double	Lz=0;

	// Read fundamental command line parameters 
		std::vector<std::string> args(argv, argv + argc);
    
		//check filename
		if (args[1] != "filename") {
			std::cout << "Arg 1 wrong";
			return 0;
		}
		else {
			filename = args[2];
			parameter_filename = args[3]; //lammps input filename
		}
		//check n_steps
		if (args[4] != "n_steps") {
			std::cout << "Arg 4 wrong";
			return 0;
		}
		else {
			n_steps = stoi(args[5]);
		}

		//check bins
		if (args[6] != "n_bins") {
			std::cout << "Arg 6 wrong";
			return 0;
		}
		else {
			n_bin = stoi(args[7]);
		}

	//Read fundamental parameters from dump file
		read_parameters(filename, n_types, n_particles, Lx, Ly, Lz);

	double	Ly_start = 0;
	double	Ly_end = Ly_start+Ly;

	//flags for the kernel
	int flag_kernel = 0;
	int flag_density = 0;
	int flag_momentum = 0; // x dir
    int flag_momentum_y = 0; // y dir
	int flag_temperature = 0;
	int flag_stressK11 = 0, flag_stressK12 = 0, flag_stressK13 = 0, flag_stressK22 = 0, flag_stressK23 = 0, flag_stressK33 = 0;
	int flag_stressV11 = 0, flag_stressV12 = 0, flag_stressV13 = 0, flag_stressV22 = 0, flag_stressV23 = 0, flag_stressV33 = 0;
	int flag_viscosity = 0;
    int flag_viscosity_bulk = 0;

	arma::Mat<double> type = arma::zeros(n_particles, n_steps);
	arma::Mat<double> x = arma::zeros(n_particles, n_steps);
	arma::Mat<double> y = arma::zeros(n_particles, n_steps);
	arma::Mat<double> z = arma::zeros(n_particles, n_steps);
	arma::Mat<double> vx = arma::zeros(n_particles, n_steps);
	arma::Mat<double> vy = arma::zeros(n_particles, n_steps);
	arma::Mat<double> vz = arma::zeros(n_particles, n_steps);
    
    arma::Mat<double> x1;
    arma::Mat<double> y1;
    arma::Mat<double> z1;
    arma::Mat<double> vx1;
    arma::Mat<double> vy1;
    arma::Mat<double> vz1;

	arma::vec bin_vec = arma::linspace<arma::vec>(Ly_start, Ly_end, n_bin);
	bin_vec.save("bin_vec.txt", arma::arma_ascii);
	arma::vec Mass = arma::ones(n_particles);
    
	double  R_cut; //characteristic length of the kernel function

	//Paramaters for lennard-jones interaction
		double Force_cut = 3.5; //cutting radius LJ
		arma::Mat<double>  epsilon = arma::ones(n_types, n_types);
		
		read_parameter_file(parameter_filename, epsilon, Force_cut);
		
    int type_wall;
    if (n_types !=1){
        type_wall = n_types - 1;
    }
    else{
        type_wall =  2;
    }
    
		double Force_cut_sq = Force_cut*Force_cut;
	///////////////////////////////////////////////////////////////

	//check kernel
	if (args[8] != "kernel") {
		std::cout << "Arg 8 wrong";
		return 0;
	}
	else {
		if (args[9] == "constant") {
			flag_kernel = 0;
		}
		else if (args[9]== "gaussian") {
			flag_kernel = 1;
		}
		else if (args[9]== "cubic") {
			flag_kernel = 2;
		}
		else {
			return 0;
		}
		R_cut = stod(args[10]);
	}

	for (size_t i_arg = 11; i_arg < args.size(); i_arg++) {
		if (args[i_arg] == "all") {
			flag_density = 1;
			flag_momentum = 1;
            flag_momentum_y = 1;
			flag_temperature = 1;
			flag_stressK11 = 1;
			flag_stressK12 = 1;
			flag_stressK13 = 1;
			flag_stressK22 = 1;
			flag_stressK23 = 1;
			flag_stressK33 = 1;
			flag_stressV11 = 1;
			flag_stressV12 = 1;
			flag_stressV13 = 1;
			flag_stressV22 = 1;
			flag_stressV23 = 1;
			flag_stressV33 = 1;
			flag_viscosity = 1;
            flag_viscosity_bulk = 1;
		}		
		if (args[i_arg] == "density") {
			flag_density = 1;
		}
		if (args[i_arg] == "momentum") {
			flag_momentum = 1;
		}
        if (args[i_arg] == "momentum_y") {
            flag_momentum_y = 1;
        }
		if (args[i_arg] == "temperature") {
            flag_momentum_y = 1;
            flag_momentum = 1;
			flag_temperature = 1;
		}
		if (args[i_arg] == "stressK11") {
			flag_stressK11 = 1;
		}
		if (args[i_arg] == "stressK12") {
            flag_momentum = 1;
			flag_stressK12 = 1;
		}
		if (args[i_arg] == "stressK13") {
			flag_stressK13 = 1;
		}
		if (args[i_arg] == "stressK22") {
			flag_stressK22 = 1;
            flag_momentum_y = 1;
		}
		if (args[i_arg] == "stressK23") {
			flag_stressK23 = 1;
		}
		if (args[i_arg] == "stressK33") {
			flag_stressK33 = 1;
		}
		if (args[i_arg] == "stressV11") {
			flag_stressV11 = 1;
		}
		if (args[i_arg] == "stressV12") {
			flag_stressV12 = 1;
		}
		if (args[i_arg] == "stressV13") {
			flag_stressV13 = 1;
		}
		if (args[i_arg] == "stressV22") {
			flag_stressV22 = 1;
		}
		if (args[i_arg] == "stressV23") {
			flag_stressV23 = 1;
		}
		if (args[i_arg] == "stressV33") {
			flag_stressV33 = 1;
		}
		if (args[i_arg] == "viscosity") {
            flag_stressV12 = 1;
            flag_stressK12 = 1;
			flag_viscosity = 1;
            flag_momentum = 1;
		}
        if (args[i_arg] == "viscosity_bulk") {
            flag_stressV22 = 1;
            flag_stressK22 = 1;
            flag_momentum_y = 1;
            flag_viscosity_bulk = 1;
        }
	}

	//read file
	std::cout << "Reading data\n";
	read_file(filename, n_steps, n_particles, type, x, y, z, vx, vy, vz);
	//std::cout << "x:\n" << x << "\n";

    arma::vec x_tmp = x.col(0);
    arma::vec type_tmp = type.col(0);
    arma::vec x1_tmp = x_tmp.elem( find(type_tmp == 0) );
    int n_particles1 = x1_tmp.n_elem;
    
    x1= reshape( x.elem( find(type == 0) ) , n_particles1, n_steps );
    y1= reshape( y.elem( find(type == 0) ) , n_particles1, n_steps );
    z1= reshape( z.elem( find(type == 0) ) , n_particles1, n_steps );
    vx1= reshape( vx.elem( find(type == 0) ) , n_particles1, n_steps );
    vy1= reshape( vy.elem( find(type == 0) ) , n_particles1, n_steps );
    vz1= reshape( vz.elem( find(type == 0) ) , n_particles1, n_steps );
    
    arma::vec Mass1 = arma::ones(x1.n_rows);
    
    //x1.print();
    
	//Compute Density
	std::cout << "Computing density \n";
	arma::vec density = arma::zeros(n_bin);
    arma::vec density1 = arma::zeros(n_bin);
	if (flag_density == 1) {
		compute_density(n_steps, n_bin, Lx, Lz, R_cut, Mass, density, y, bin_vec, flag_kernel);
		density.save("density.txt", arma::arma_ascii);
        
        compute_density(n_steps, n_bin, Lx, Lz, R_cut, Mass1, density1, y1, bin_vec, flag_kernel);
        density1.save("density1.txt", arma::arma_ascii);
	}

	//Compute Velocity in x direction
	arma::vec momentum = arma::zeros(n_bin);
	arma::vec velocity;
    
    arma::vec momentum1 = arma::zeros(n_bin);
    arma::vec velocity1;
    
	arma::vec zero_velocity = arma::zeros(n_bin);
	if (flag_momentum == 1) {
		std::cout << "Computing momentum and velocity in x direction\n";
		compute_momentum(n_steps, n_bin, Lx, Lz, R_cut, vx, momentum, y, bin_vec, flag_kernel);
		velocity = momentum / density;
		momentum.save("momentum.txt", arma::arma_ascii);
		velocity.save("velocity.txt", arma::arma_ascii);
        
        compute_momentum(n_steps, n_bin, Lx, Lz, R_cut, vx1, momentum1, y1, bin_vec, flag_kernel);
        velocity1 = momentum1 / density1;
        momentum1.save("momentum1.txt", arma::arma_ascii);
        velocity1.save("velocity1.txt", arma::arma_ascii);
        
	}
    
    //Compute Velocity in y direction
    arma::vec momentum_y = arma::zeros(n_bin);
    arma::vec velocity_y;
    if (flag_momentum_y == 1) {
        std::cout << "Computing momentum and velocity in y direction\n";
        compute_momentum_y(n_steps, n_bin, Lx, Lz, R_cut, vy, momentum_y, y, bin_vec, flag_kernel);
        velocity_y = momentum_y / density;
        momentum_y.save("momentum_y.txt", arma::arma_ascii);
        velocity_y.save("velocity_y.txt", arma::arma_ascii);
    }

	//Compute Temperature
	if (flag_temperature == 1) {
		std::cout << "Computing temperature \n";
		arma::vec temp_tmp = arma::zeros(n_bin);
		compute_temperature(n_steps, n_bin, Lx, Lz, R_cut, vx, vy, vz, velocity,velocity_y, temp_tmp, y, bin_vec, flag_kernel);
		arma::vec temperature = temp_tmp / density;
		temperature.save("temperature.txt", arma::arma_ascii);
	}

	//Compute Kinetic Stresses
	if (flag_stressK11 == 1) {
		std::cout << "Computing kinetic stress 11 \n";
		arma::vec stressK = arma::zeros(n_bin);
        arma::vec stressK1 = arma::zeros(n_bin);
		compute_stressK(n_steps, n_bin, Lx, Lz, R_cut, vx, vx, velocity, velocity, stressK, y, bin_vec, flag_kernel);
        compute_stressK(n_steps, n_bin, Lx, Lz, R_cut, vx1, vx1, velocity1, velocity1, stressK1, y1, bin_vec, flag_kernel);
		stressK.save("stressK11.txt", arma::arma_ascii);
        stressK1.save("stressK11_1.txt", arma::arma_ascii);
	}
	
	arma::vec stressK12 = arma::zeros(n_bin);
    arma::vec stressK12_1 = arma::zeros(n_bin);
	if (flag_stressK12 == 1) {
		std::cout << "Computing kinetic stress 12 \n";
		compute_stressK(n_steps, n_bin, Lx, Lz, R_cut, vx, vy, velocity, zero_velocity, stressK12, y, bin_vec, flag_kernel);
        compute_stressK(n_steps, n_bin, Lx, Lz, R_cut, vx1, vy1, velocity1, zero_velocity, stressK12_1, y1, bin_vec, flag_kernel);
		stressK12.save("stressK12.txt", arma::arma_ascii);
        stressK12_1.save("stressK12_1.txt", arma::arma_ascii);
	}
	if (flag_stressK13 == 1) {
		std::cout << "Computing kinetic stress 13 \n";
		arma::vec stressK = arma::zeros(n_bin);
        arma::vec stressK1 = arma::zeros(n_bin);
		compute_stressK(n_steps, n_bin, Lx, Lz, R_cut, vx, vz, velocity, zero_velocity, stressK, y, bin_vec, flag_kernel);
        compute_stressK(n_steps, n_bin, Lx, Lz, R_cut, vx1, vz1, velocity1, zero_velocity, stressK1, y1, bin_vec, flag_kernel);
        stressK.save("stressK13.txt", arma::arma_ascii);
		stressK1.save("stressK13_1.txt", arma::arma_ascii);
	}
    

	if (flag_stressK22 == 1) {
		std::cout << "Computing kinetic stress 22 \n";
		arma::vec stressK = arma::zeros(n_bin);
        arma::vec stressK1 = arma::zeros(n_bin);
		compute_stressK(n_steps, n_bin, Lx, Lz, R_cut, vy, vy, velocity_y, velocity_y, stressK, y, bin_vec, flag_kernel);
        compute_stressK(n_steps, n_bin, Lx, Lz, R_cut, vy1, vy1, velocity_y, velocity_y, stressK1, y1, bin_vec, flag_kernel);
		stressK.save("stressK22.txt", arma::arma_ascii);
        stressK1.save("stressK22_1.txt", arma::arma_ascii);
	}
	if (flag_stressK23 == 1) {
		std::cout << "Computing kinetic stress 23 \n";
		arma::vec stressK = arma::zeros(n_bin);
        arma::vec stressK1 = arma::zeros(n_bin);
		compute_stressK(n_steps, n_bin, Lx, Lz, R_cut, vy, vz, zero_velocity, zero_velocity, stressK, y, bin_vec, flag_kernel);
        compute_stressK(n_steps, n_bin, Lx, Lz, R_cut, vy1, vz1, zero_velocity, zero_velocity, stressK1, y1, bin_vec, flag_kernel);
        stressK.save("stressK23.txt", arma::arma_ascii);
		stressK1.save("stressK23_1.txt", arma::arma_ascii);
	}
	if (flag_stressK33 == 1) {
		std::cout << "Computing kinetic stress 33 \n";
		arma::vec stressK = arma::zeros(n_bin);
        arma::vec stressK1 = arma::zeros(n_bin);
		compute_stressK(n_steps, n_bin, Lx, Lz, R_cut, vz, vz, zero_velocity, zero_velocity, stressK, y, bin_vec, flag_kernel);
        compute_stressK(n_steps, n_bin, Lx, Lz, R_cut, vz1, vz1, zero_velocity, zero_velocity, stressK1, y1, bin_vec, flag_kernel);
        stressK.save("stressK33.txt", arma::arma_ascii);
		stressK1.save("stressK33_1.txt", arma::arma_ascii);
	}


	//Compute Potential Stresses
	arma::vec stressV12 = arma::zeros(n_bin);
    arma::vec stressV12_1 = arma::zeros(n_bin);
    
	if (flag_stressV11 == 1 || flag_stressV12 == 1 || flag_stressV13 == 1 || flag_stressV22 == 1 || flag_stressV23 == 1 || flag_stressV33 == 1) {
		std::cout << "Computing potential stress tensor \n";

		arma::vec stressV11 = arma::zeros(n_bin);
		arma::vec stressV13 = arma::zeros(n_bin);
		arma::vec stressV22 = arma::zeros(n_bin);
		arma::vec stressV23 = arma::zeros(n_bin);
		arma::vec stressV33 = arma::zeros(n_bin);
        
        arma::vec stressV11_1 = arma::zeros(n_bin);
        arma::vec stressV13_1 = arma::zeros(n_bin);
        arma::vec stressV22_1 = arma::zeros(n_bin);
        arma::vec stressV23_1 = arma::zeros(n_bin);
        arma::vec stressV33_1 = arma::zeros(n_bin);
        
		//double stress_virial = 0.0;
		double force_wall_h = 0.0;
		double force_wall_v = 0.0;
        double erf_factor = 1 / (std::sqrt(2.0) * R_cut);
		
        for (size_t i_step = 0; i_step < n_steps; i_step++) {
			double pair_energy = 0.0;
			

			for (size_t i = 1; i < n_particles; i++) {
				double x_i = x(i, i_step);
				double y_i = y(i, i_step);
				double z_i = z(i, i_step);

				for (size_t j = 0; j < n_particles ; j++) {
					double y_j = y(j, i_step);
					//Minimim image convention due to periodic BC in x and z directions
					
					double distx = x(j, i_step) - x_i;
					if (distx > 0.5*Lx ) {
						distx = distx - Lx;
					}
					else if (distx < -0.5*Lx ) {
						distx = distx + Lx;
					}

                    double disty = y_j - y_i;
//                    if (disty > 0.5*Ly) {
//                        disty = disty -  Ly;
//                        y_j = y_j - Ly;
//                    }
//                    else if (disty < -0.5 * Ly) {
//                        disty = disty + Ly;
//                        y_j= y_j + Ly;
//                    }

					double distz = z(j, i_step) - z_i;
					if (distz > 0.5*Lz) {
						distz = distz - Lz;
					}
					else if (distz < -0.5*Lz) {
						distz = distz + Lz;
					}

					double r_sq = distx*distx + disty*disty + distz*distz;
                    //r_sq > 0.0 &&
                    
					if ( i!=j && r_sq > 0 && r_sq <= Force_cut_sq &&  type(i) ==0  ) {
						double r = std::sqrt(r_sq);
						double r6i = 1.0 / (r_sq*r_sq*r_sq);
						double force = epsilon(type(i,i_step), type(j,i_step)) *r6i * (48.0 * r6i - 24.0) / r; 

						if (type(i) == type_wall && type(j) != type_wall) {
							force_wall_h += force*distx / r / (Lx * Lz) /n_steps ;
							force_wall_v += force*disty / r / (Lx * Lz) / n_steps;
						}
						else if (type(j) == type_wall && type(i) != type_wall) {
							force_wall_h -= force*distx / r / (Lx * Lz)  /n_steps;
							force_wall_v -= force*disty / r / (Lx * Lz) / n_steps;
						}

						//stress_virial -= (1.0 / (3.0* Lx * Ly * Lz)) * force * r / n_steps;
						//pair_energy += 4.0 * r6i * (r6i - 1.0) / n_particles ;
                            
						for (int i_bin = 0; i_bin < n_bin; i_bin++) {
							double y_bin = bin_vec(i_bin);
							if (y_i > y_bin - R_cut - Force_cut && y_j < y_bin + R_cut + Force_cut && y_i> y_bin - R_cut - Force_cut && y_i < y_bin + R_cut + Force_cut) {
								//compute bond function
								double bond_function;
								if (disty == 0 )  {
                                    //printf("P: %d %d", (int) type(i) , (int) type(j) );
									bond_function = 0;
                                }
								else if (flag_kernel == 0) {
									bond_function = - 1 / (2.0 * R_cut * abs(disty)) * std::max(0.0, std::min ( y_bin + R_cut, (std::max(y_i, y_j))) - std::max(y_bin - R_cut, std::min(y_i, y_j)) );
								}
								else if (flag_kernel == 1) {
									bond_function = 1 / (2.0 * disty ) *(std::erf(erf_factor*(y_bin - y_j)) - std::erf(erf_factor*(y_bin - y_i)));
								}
                                else if (flag_kernel == 2) {
                                    bond_function = 1 / ( disty * R_cut) * ( cubic_bond_function( y_bin - y_j ,R_cut ) - cubic_bond_function( y_bin - y_i , R_cut ) ) ;
                                }

								if (flag_stressV11 == 1) {
									stressV11(i_bin) += 0.5 / (Lx * Lz) * distx * (force * distx / r) * bond_function / n_steps;
                                    stressV11_1(i_bin) += 0.5 / (Lx * Lz) * distx * (force * distx / r) * bond_function / n_steps;
								}
								if (flag_stressV12 == 1) {
									stressV12(i_bin) += 0.5 / (Lx * Lz) * distx * (force * disty / r) * bond_function / n_steps;
                                    stressV12_1(i_bin) += 0.5 / (Lx * Lz) * distx * (force * disty / r) * bond_function / n_steps;
								}
								if (flag_stressV13 == 1) {
									stressV13(i_bin) += 0.5 / (Lx * Lz) * distx * (force * distz / r) * bond_function / n_steps;
                                    stressV13_1(i_bin) += 0.5 / (Lx * Lz) * distx * (force * distz / r) * bond_function / n_steps;
								}
								if (flag_stressV22 == 1) {
									stressV22(i_bin) += 0.5 / (Lx * Lz) * disty * (force * disty / r) * bond_function / n_steps;
                                    stressV22_1(i_bin) += 0.5 / (Lx * Lz) * disty * (force * disty / r) * bond_function / n_steps;
								}
								if (flag_stressV23 == 1) {
									stressV23(i_bin) += 0.5 / (Lx * Lz) * disty * (force * distz / r) * bond_function / n_steps;
                                    stressV23_1(i_bin) += 0.5 / (Lx * Lz) * disty * (force * distz / r) * bond_function / n_steps;
								}
								if (flag_stressV33 == 1) {
									stressV33(i_bin) += 0.5 / (Lx * Lz) * distz * (force * distz / r) * bond_function / n_steps;
                                    stressV33_1(i_bin) += 0.5 / (Lx * Lz) * distz * (force * distz / r) * bond_function / n_steps;
								}
							}//end if
						}//end for
                    } else if ( i!=j && r_sq > 0 && r_sq <= Force_cut_sq )
                        {
                        double r = std::sqrt(r_sq);
                        double r6i = 1.0 / (r_sq*r_sq*r_sq);
                        double force = epsilon(type(i,i_step), type(j,i_step)) *r6i * (48.0 * r6i - 24.0) / r;
                        
                        if (type(i) == type_wall && type(j) != type_wall) {
                            force_wall_h += force*distx / r / (Lx * Lz) /n_steps ;
                            force_wall_v += force*disty / r / (Lx * Lz) / n_steps;
                        }
                        else if (type(j) == type_wall && type(i) != type_wall) {
                            force_wall_h -= force*distx / r / (Lx * Lz)  /n_steps;
                            force_wall_v -= force*disty / r / (Lx * Lz) / n_steps;
                        }
                        
                        //stress_virial -= (1.0 / (3.0* Lx * Ly * Lz)) * force * r / n_steps;
                        //pair_energy += 4.0 * r6i * (r6i - 1.0) / n_particles ;
                        
                        for (int i_bin = 0; i_bin < n_bin; i_bin++) {
                            double y_bin = bin_vec(i_bin);
                            if (y_i > y_bin - R_cut - Force_cut && y_j < y_bin + R_cut + Force_cut && y_i> y_bin - R_cut - Force_cut && y_i < y_bin + R_cut + Force_cut) {
                                //compute bond function
                                double bond_function;
                                if (disty == 0 )  {
                                    //printf("P: %d %d", (int) type(i) , (int) type(j) );
                                    bond_function = 0;
                                }
                                else if (flag_kernel == 0) {
                                    bond_function = - 1 / (2.0 * R_cut * abs(disty)) * std::max(0.0, std::min ( y_bin + R_cut, (std::max(y_i, y_j))) - std::max(y_bin - R_cut, std::min(y_i, y_j)) );
                                }
                                else if (flag_kernel == 1) {
                                    bond_function = 1 / (2.0 * disty ) *(std::erf(erf_factor*(y_bin - y_j)) - std::erf(erf_factor*(y_bin - y_i)));
                                }
                                else if (flag_kernel == 2) {
                                    bond_function = 1 / ( disty * R_cut) * ( cubic_bond_function( y_bin - y_j ,R_cut ) - cubic_bond_function( y_bin - y_i , R_cut ) ) ;
                                }
                                
                                if (flag_stressV11 == 1) {
                                    stressV11(i_bin) += 0.5 / (Lx * Lz) * distx * (force * distx / r) * bond_function / n_steps;
                                }
                                if (flag_stressV12 == 1) {
                                    stressV12(i_bin) += 0.5 / (Lx * Lz) * distx * (force * disty / r) * bond_function / n_steps;
                                }
                                if (flag_stressV13 == 1) {
                                    stressV13(i_bin) += 0.5 / (Lx * Lz) * distx * (force * distz / r) * bond_function / n_steps;
                                }
                                if (flag_stressV22 == 1) {
                                    stressV22(i_bin) += 0.5 / (Lx * Lz) * disty * (force * disty / r) * bond_function / n_steps;
                                }
                                if (flag_stressV23 == 1) {
                                    stressV23(i_bin) += 0.5 / (Lx * Lz) * disty * (force * distz / r) * bond_function / n_steps;
                                }
                                if (flag_stressV33 == 1) {
                                    stressV33(i_bin) += 0.5 / (Lx * Lz) * distz * (force * distz / r) * bond_function / n_steps;
                                }
                            }//end if
                        }//end for
                    }
                    
				} //end for j
			}//end for i
				//std::cout << "Average Pair Energy/particle" << pair_energy << "\n";
				//std::cout << "Stress Virial" << stress_virial*n_steps / (i_step + 1) << "\n";
				//std::cout << "Force wall" << force_wall*n_steps / (i_step + 1) << "\n";
			} //end for nsteps
			//std::cout << stress_virial << "\n";

			if (flag_stressV11 == 1) {
				stressV11.save("stressV11.txt", arma::arma_ascii);
                stressV11_1.save("stressV11_1.txt", arma::arma_ascii);
			}
			if (flag_stressV12 == 1) {
				stressV12.save("stressV12.txt", arma::arma_ascii);
                stressV12_1.save("stressV12_1.txt", arma::arma_ascii);
			}
			if (flag_stressV13 == 1) {
				stressV13.save("stressV13.txt", arma::arma_ascii);
                stressV13_1.save("stressV13_1.txt", arma::arma_ascii);
			}
			if (flag_stressV22 == 1) {
				stressV22.save("stressV22.txt", arma::arma_ascii);
                stressV22_1.save("stressV22_1.txt", arma::arma_ascii);
			}
			if (flag_stressV23 == 1) {
				stressV23.save("stressV23.txt", arma::arma_ascii);
                stressV23_1.save("stressV23_1.txt", arma::arma_ascii);
			}
			if (flag_stressV33 == 1) {
				stressV33.save("stressV33.txt", arma::arma_ascii);
                stressV33_1.save("stressV33_1.txt", arma::arma_ascii);
			}

			//write Stress at wall on file
			std::ofstream myfile("WallStress.txt");
			if (myfile.is_open())
			{
				myfile << "Horizontal stress  " << force_wall_h << "\n";
				myfile << "Vertical stress  " << force_wall_v << "\n";
				myfile.close();
			}

		} //end if

		//Compute viscosity
		if (flag_viscosity == 1 && flag_momentum == 1) {
			std::cout << "Computing viscosity \n";
			arma::vec velocity_grad = arma::zeros(n_bin);
			arma::vec viscosity = arma::zeros(n_bin);
            arma::vec viscosity1 = arma::zeros(n_bin);
            arma::vec density_derivative = arma::zeros(n_bin);
            arma::vec momentum_derivative = arma::zeros(n_bin);
            
//            if (flag_kernel == 0 || flag_kernel == 2) {
            if (flag_kernel == 0) {
                for (size_t i_bin = 1; i_bin < n_bin - 1; i_bin++) {
                    velocity_grad(i_bin) = ((momentum(i_bin + 1) - momentum(i_bin - 1)) * density(i_bin) - (density(i_bin + 1) - density(i_bin - 1)) * momentum(i_bin)) / ((bin_vec(i_bin + 1)-bin_vec(i_bin - 1))* density(i_bin)*density(i_bin));
                    if (velocity_grad(i_bin) != 0.0) {
                        viscosity(i_bin) = (stressK12(i_bin) + stressV12(i_bin)) / velocity_grad(i_bin);
                        viscosity1(i_bin) = (stressK12_1(i_bin) + stressV12_1(i_bin)) / velocity_grad(i_bin);
                    }
                    
                    density_derivative(i_bin) = (density(i_bin + 1) - density(i_bin - 1))/ (bin_vec(i_bin + 1)-bin_vec(i_bin - 1));
                    momentum_derivative(i_bin) = (momentum(i_bin + 1) - momentum(i_bin - 1))/ (bin_vec(i_bin + 1)-bin_vec(i_bin - 1));
                }
                density_derivative.save("density_derivative.txt", arma::arma_ascii);
                momentum_derivative.save("momentum_derivative.txt", arma::arma_ascii);
            }else{
                compute_density_derivative(n_steps, n_bin, Lx, Lz, R_cut, Mass, density_derivative, y, bin_vec, flag_kernel);
                
                compute_momentum_derivative(n_steps, n_bin, Lx, Lz, R_cut, vx, momentum_derivative, y, bin_vec, flag_kernel);
                
                velocity_grad = (momentum_derivative - velocity % density_derivative) / density;
                
                viscosity= (stressK12+stressV12) / velocity_grad;
                viscosity1= (stressK12_1+stressV12_1) / velocity_grad;
                
                density_derivative.save("density_derivative.txt", arma::arma_ascii);
                momentum_derivative.save("momentum_derivative.txt", arma::arma_ascii);
            }
            
			viscosity.save("viscosity.txt", arma::arma_ascii);
            viscosity1.save("viscosity1.txt", arma::arma_ascii);
			velocity_grad.save("velocity_grad.txt", arma::arma_ascii);
		}
    
		
	printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
		//std::cin.get();
	return 0;
}
