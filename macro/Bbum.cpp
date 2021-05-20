#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <ctime>
#include <random>
#include <list>
#include <thread>
#include "vecFunc.cpp"
//#include <Python.h>
//using std::vector;

//g++ -pthread -o Bbum.out Bbum.cpp


class Wall {
private:
	double wall_cell_size = 0.5;
	long int wall_number_of_particles = 150;
	long int wall_number_of_particles2 = std::pow(wall_number_of_particles, 2);
	//std::vector<double> v(2);
	std::vector<std::vector<double>> data_wall = std::vector<std::vector<double>> (wall_number_of_particles2, std::vector<double> (2));
	//std::vector<double> data_wall (wall_number_of_particles2);
	double wall_coord; // = - 2*R_intr
	long int wall_rude_coord;// = wall_coord//R_intr + field_half_size


	void generate_wall(){
		for (int i = 0; i<wall_number_of_particles; i++){
			for (int j = 0; j<wall_number_of_particles; j++){
				std::cout << i << "\t" << j << "\n";
				data_wall[wall_number_of_particles*i+j][0] = (i - wall_number_of_particles/2)*wall_cell_size;
				data_wall[wall_number_of_particles*i+j][1] = (j - wall_number_of_particles/2)*wall_cell_size;
			}
		}

	}


public:
	Wall(double wc){

		//data_wall = std::vector<std::vector<double>> (wall_number_of_particles2, std::vector<double> (2));
		generate_wall();
		wall_coord = wc;
		output();
	}

	void output(){
		std::cout << "Введите имя стенофайла:\n";
		std::string name;
		std::cin >> name;
		std::ofstream file("results/wall" + name + ".xyz");
		file << wall_number_of_particles2 << "\n\n";
		for (int j=0; j<wall_number_of_particles2; j++){
			file << "10\t" << wall_coord << "\t" << data_wall[j][0] << "\t" << data_wall[j][1] << "\n";
		}

	}


};



class Calc {
private:
	long int number_of_particles;

	float R_intr = 3; // примерно 7 равновесных радиусов 50 \approx 7^2
	float R_intr2 = std::pow(R_intr, 2);
	double tau = 0.0005;
	double tau2 = std::pow(tau, 2);
	long int max_step = 2000;
	double delay_frame = 0.002;
	long int field_half_size = (long int) (100/R_intr);
	long int field_size = 2*field_half_size+1;
	double force_size = 10;

	long int step = 1;

	double wall_coord = - 4*R_intr;
	long int wall_rude_coord = ((long int) wall_coord/R_intr) + field_half_size;

	//std::vector<std::vector<std::vector<double>>> data;
	std::vector<std::vector<long int>> rude_coords;
	std::vector<std::vector<std::vector<std::list<long int>>>> field;

	std::vector<std::vector<double>> * pp_coords;
	std::vector<std::vector<double>> * pc_coords;

	std::vector<std::vector<double>> past_coords;
	std::vector<std::vector<double>> cur_coords;



	time_t time = std::time(0);

	double percent = 0;

	int number_of_cores = std::thread::hardware_concurrency();

	void generate_coords(){
		rude_coords = std::vector<std::vector<long int>> (number_of_particles, std::vector<long int> (3));

		long int size = (long int) round(std::pow(number_of_particles, (1.0/3.0)));

		//coords = np.zeros((number_of_particles, 3))
		for (long int i=0; i < size; i++){
			for (long int j=0; j < size; j++){
				for (long int k=0; k < size; k++){
					(*pp_coords)[size*(size*i+j)+k] = {(double)i, (double) j-size/2, (double) k-size/2};
				}
			}
		}
		//coords = 0.9*coords
		//coords += 0.1*np.random.rand(number_of_particles, 3) - 0.05

		std::vector<std::vector<double>> velocities = std::vector<std::vector<double>> (number_of_particles, std::vector<double> (3));
		
		for (int i=0; i<number_of_particles; i++){
			random_v(velocities[i]);
			velocities[i] += -0.5;
			//std::cout << velocities[i][0] << ' ' << velocities[i][1] << ' ' << velocities[i][2] << '\n';
		}
		
		if (wall_status){
			for (int i=0; i<number_of_particles; i++){
				velocities[i][0] += (wall_coord-size/2)*2.5;///(tau*10000)*2.5; //max_step здесь принято за 1E5
			}
		}
		velocities = 2*velocities;

		(*pc_coords) = (*pp_coords) + tau*velocities;
		//std::cout << data[1][0][0] - data[0][0][0] << '\n';
	}

	
	void generate_field(){

		field.clear();
		field.resize(field_size);
		for (int i=0; i<field_size; i++){
			field[i].resize(field_size);
			for (int j=0; j<field_size; j++){
				field[i][j].resize(field_size);
			}

		}
		rude_coords.resize(number_of_particles);
		rude_coords =  rounder((*pc_coords)/R_intr + field_half_size); //т.к rude_coord -- целочисленный, то деление автоматически тоже целочисленное
		for (int i=0; i<number_of_particles; i++){
			if ((rude_coords[i][0] >= 0) and (rude_coords[i][1] >= 0) and (rude_coords[i][2] >= 0) 
			and (rude_coords[i][0] < field_size) and (rude_coords[i][1] < field_size) and (rude_coords[i][2] < field_size)){
				field[rude_coords[i][0]][rude_coords[i][1]][rude_coords[i][2]].push_back(i);
			}
		}
	}

	std::vector<double> force(std::vector<double> vec){
		double R = sum(power(vec,2)); //квадрат расстояния
		if (R > R_intr2)
			return std::vector<double> (3);
		double f =  force_size*(std::pow(R, -7) - std::pow(R,-4)); // \sigma = (1/2)^(1/6)
		return f*vec;
	}

	double reflection(long int number){
		double r = std::abs((*pc_coords)[number][0]-wall_coord);
		return std::pow(r, -13);
	}

	/*
	std::vector<std::vector<double>> velocities(long int stp){
		if (stp == 0){
			stp = 1;
		}
		return (data[stp] - data[stp-1])/tau; 
	}
	*/


	std::vector<std::vector<double>> accelerations(){
		//std::vector<>coords = (*pc_coords)
		std::vector<std::vector<double>> accelerations = std::vector<std::vector<double>>(number_of_particles, std::vector<double>(3));
		generate_field();
		for (long int i =0; i< number_of_particles; i++){
			if (wall_status and (std::abs(rude_coords[i][0] - wall_rude_coord) < 1)){
				accelerations[i][0] += reflection(i);
			}
			if ((rude_coords[i][0] >= 1) and (rude_coords[i][1] >= 1) and (rude_coords[i][2] >= 1) 
			and (rude_coords[i][0] < field_size-1) and (rude_coords[i][1] < field_size-1) and (rude_coords[i][2] < field_size-1)){
				//if (field[rude_coords[i][0]][rude_coords[i][1]][rude_coords[i][2]].front() != i){
					//std::cout << "SOS!" << field[rude_coords[i][0]][rude_coords[i][1]][rude_coords[i][2]].front() << '\n';
				//}
				field[rude_coords[i][0]][rude_coords[i][1]][rude_coords[i][2]].pop_front();
				for (long int x = -1; x<2; x++){
					for (long int y = -1; y<2; y++){
						for (long int z = -1; z<2; z++){
							std::list <long int> :: iterator cell;
							for (cell = field[rude_coords[i][0]+x][rude_coords[i][1]+y][rude_coords[i][2]+z].begin();
							    cell != field[rude_coords[i][0]+x][rude_coords[i][1]+y][rude_coords[i][2]+z].end(); cell++){
								//print("F",i,j)
								std::vector<double> forceic = force(((*pc_coords)[i]-(*pc_coords)[*cell]));
								accelerations[i] += forceic;
								accelerations[*cell] -= forceic;
							}
						}
					}
				}
			}
		}
		//print("Acc:\n")
		//print(accelerations)
		return accelerations;
	}

	std::vector<std::vector<double>> m_accelerations(){
		std::vector<std::vector<double>> accelerations = std::vector<std::vector<double>>(number_of_particles, std::vector<double>(3));
		std::vector<std::vector<std::vector<double>>> local_accelerations = std::vector<std::vector<std::vector<double>>>(number_of_cores, std::vector<std::vector<double>>(number_of_particles, std::vector<double>(3)));
		generate_field();
		std::thread* helper = new std::thread[number_of_cores];
    	for (int j=0; j<number_of_cores; j++){
        	helper[j] = std::thread(&Calc::m_calc_acc, this, j, std::ref(local_accelerations[j]));
    	}
    	for (int i=0; i<number_of_cores; i++){
        	helper[i].join();
   		}
   		for (int i=0; i<number_of_cores; i++){
   			accelerations += local_accelerations[i];
   		}
   		delete [] helper;
   		return accelerations;
   	}

	
	void m_calc_acc(long int beg, std::vector<std::vector<double>> & accelerations){
		//std::vector<>coords = (*pc_coords)
		//long int NP = number_of_particles, NC = number_of_cores;
		for (long int i=beg; i<number_of_particles; i+=number_of_cores){ //возможно тут стоит создать новую переменную, чтобы все подряд не обращались к number_of_particles
			if (wall_status and (std::abs(rude_coords[i][0] - wall_rude_coord) < 1)){
				accelerations[i][0] += reflection(i);
				//acc[0] += reflection(i);
			}
			if ((rude_coords[i][0] >= 1) and (rude_coords[i][1] >= 1) and (rude_coords[i][2] >= 1) 
			and (rude_coords[i][0] < field_size-1) and (rude_coords[i][1] < field_size-1) and (rude_coords[i][2] < field_size-1)){
				//if (field[rude_coords[i][0]][rude_coords[i][1]][rude_coords[i][2]].front() != i){
					//std::cout << "SOS!" << field[rude_coords[i][0]][rude_coords[i][1]][rude_coords[i][2]].front() << '\n';
				//}
				for (long int x = -1; x<2; x++){
					for (long int y = -1; y<2; y++){
						for (long int z = -1; z<2; z++){
							std::list <long int> :: iterator cell;
							std::list <long int> :: iterator cell_end = field[rude_coords[i][0]+x][rude_coords[i][1]+y][rude_coords[i][2]+z].end();
							for (cell = field[rude_coords[i][0]+x][rude_coords[i][1]+y][rude_coords[i][2]+z].begin();
							    cell != cell_end; cell++){
								//print("F",i,j)
								if (i < *cell){
									std::vector<double> acc;
									acc.resize(3);
									acc = force(((*pc_coords)[i]-(*pc_coords)[*cell]));
									accelerations[i] += acc;
									accelerations[*cell] -= acc;
								}
							}
						}
					}
				}
			}
		}
	}

	void evolve(std::ofstream & file, std::ofstream &calcs){
		//немного медленнее, но оперативка O(N)


		step = 1;
		percent = 0;
		time = std::time(0);

		past_coords = std::vector<std::vector<double>> (number_of_particles, std::vector<double> (3));
		cur_coords = std::vector<std::vector<double>> (number_of_particles, std::vector<double> (3));
		pp_coords = & past_coords;
		pc_coords = & cur_coords;
		std::vector<std::vector<double>> next_coords(number_of_particles, std::vector<double>(3));
		double T = delay_frame;

		generate_coords();

		double vel_of_CM_0;
		if (target == "VelCM"){
			std::cout << "N = " << number_of_particles << "\n";
			vel_of_CM_0 = velocity_of_CM();
		}
		while(step < max_step-1){

			if (multicore_status){
				next_coords = 2*cur_coords - past_coords + tau2*m_accelerations();
			}
			else {
				next_coords = 2*cur_coords - past_coords + tau2*accelerations();
			}

			past_coords = cur_coords;
			cur_coords = next_coords;
			step ++;
			loading();
			/*
			double pr = step*1000/max_step/10.0;
			if (percent < pr){
				percent = pr;
				load << pr << "%\t" << std::time(0) - time << "s\n";
			}
			*/

			T += tau;
			if (T > delay_frame){
				if (saving_status){
					file << number_of_particles << "\n\n";
					//std::cout << "Peace\n";
					for (long int j=0; j<number_of_particles; j++){
						file << "1\t" << cur_coords[j][0] << "\t" << cur_coords[j][1] << "\t" << cur_coords[j][2] << "\n";
					}
				}
				if (target == "Energy"){
					calcs << energy() << '\n';
				}
				if (target == "Maxwell"){
					calcs << number_of_particles << "\n";
					for (long int j=0; j<number_of_particles; j++){
						calcs << sum(power(cur_coords[j] - past_coords[j], 2))/tau2 << "\n";
					}
				}
				T = 0;
			}

		}
		if (target == "VelCM"){
			calcs << velocity_of_CM()/vel_of_CM_0 << '\n';
		}

	}


public:
	bool wall_status = false;
	bool multicore_status = true;
	bool saving_status = true;
	std::string target = "";

	Calc(int N){
		
		number_of_particles = std::pow((long int)(std::pow(N+1,(1.0/3.0))), 3);

		
	}

	/*
	void evolve(){
		//быстрый, но оперативка расходуется как O(S*N)

		if (((long long) number_of_particles)*((long long) max_step) < 5E8){

			step = 1;
			percent = 0;
			time = std::time(0);

			data = std::vector<std::vector<std::vector<double>>>(max_step, std::vector<std::vector<double>> (number_of_particles, std::vector<double> (3)));
			pp_coords = & (data[0]);
			pc_coords = & (data[1]);

			generate_coords();

			while(step < max_step-1){
				if (multicore_status){
					data[step+1] = 2*data[step] - data[step-1] + tau2*m_accelerations();
				}
				else {
					data[step+1] = 2*data[step] - data[step-1] + tau2*accelerations();
				}
				step ++;
				loading();
				pc_coords = & data[step];
				pp_coords = & data[step-1];
			}
		}
		else {
			std::cout << "Используйте slow_evolve()\n";
		}
		//output()
		//output_wall()

		//check_energies()
		
		//print(field)

		//E = energy(max_step-1)[2]
		//E_0 = energy(2)[2]
		//print(f"E = {E}, E_0 = {E_0}, Eps = {E/E_0 - 1}")

		//std::cout << "multicores = " << multicore_status << ' ' << std::time(0) - time_evolve << "s \n";
		
	}
	*/
	void calculus(std::string tar){
		target = tar;
		std::cout << "Введите имя файла:\n";
		std::string name; //= "";
		std::cin >> name;
		std::ofstream file;
		std::ofstream calcs;
		saving_status = true;
		if (saving_status){
			file = std::ofstream("results/bum" + name + ".xyz");
		}
		if (target != ""){
			calcs = std::ofstream("results/" + target + name + ".xyz");
		}
		if ((not saving_status) and (target == "")){
			std::cout << "Ничего не записывается\n";
		}
		if (target == "VelCM"){
			saving_status = false;
			for (int i = 1; i <15; i++){
				number_of_particles = std::pow(i, 3);
				evolve(file, calcs);
			}
			saving_status = true;
			number_of_particles = std::pow(15, 3);
			evolve(file, calcs);
		}
		else {
			evolve(file, calcs);
		}
		if (saving_status){
			file.close();
		}
		if (target != ""){
			calcs.close();
		}


		target = "";

	}
	
	/*

	void output(){
		std::cout << "Введите имя файла:\n";
		std::string name;
		std::cin >> name;
		std::ofstream file("results/bum" + name + ".xyz");
		double T = delay_frame;
		for (long int i=0; i<max_step; i++){
			T += tau;
			//std::cout << T << '\n';
			if (T >= delay_frame){
				file << number_of_particles << "\n\n";
				//std::cout << "Peace\n";
				for (long int j=0; j<number_of_particles; j++){
					file << "1\t" << data[i][j][0] << "\t" << data[i][j][1] << "\t" << data[i][j][2] << "\n";
				}
				T = 0;
			} 
		}
		file.close();
	}
	*/


	double velocity_of_CM(){
		//std::vector<double> velocity_of_CM(3);
		double vel_of_CM = 0;
		for (long int i=0; i<number_of_particles; i++){
			vel_of_CM += cur_coords[i][0] - past_coords[i][0];
		}
		return vel_of_CM/tau;
	}
	/*

	double end_analysis(){
		std::cout << "Потрачено времени : " << std::time(0) - time << "s\n";
		//double E_beg = energy(0), E_end = energy(max_step-1);  
		//std::cout << "E beg = " << E_beg << "\tE end = " << E_end << "\teps = " << (E_beg-E_end)/std::abs(E_beg) << "\n";
		double v_beg = velocity_of_CM(0)[0], v_end = velocity_of_CM(max_step-1)[0]; 
		std::cout << "v beg = " << v_beg << "\tv end = " << v_end << "\teps = " << (v_beg+v_end)/std::abs(v_beg) << "\n";
	}
	*/
	void output_wall(){
		if (wall_status){
			Wall w = Wall(wall_coord);
		}
		else {
			std::cout << "Стену нужно возвести\n";
		}
	}

	double energy(){
		//print("slow call")
		double U = 0;

		std::vector<double> local_energies(number_of_cores);
		std::thread* helper = new std::thread[number_of_cores];
    	for (int j=0; j<number_of_cores; j++){
        	helper[j] = std::thread(&Calc::m_calc_energy, this, j, std::ref(local_energies[j]));
    	}
    	for (int i=0; i<number_of_cores; i++){
        	helper[i].join();
   		}
   		delete [] helper;
   		U  = force_size*sum(local_energies)/6;


		double T = 1/2*sum(power(cur_coords-past_coords, 2))/tau2;
		return U+T;
	}

	void m_calc_energy(int beg, double & local_energies){
		//std::vector<>coords = (*pc_coords)
		//long int NP = number_of_particles, NC = number_of_cores;
		double U_piece = (1/(2*std::pow(R_intr, 6))-1)/std::pow(R_intr, 6);
		for (long int i=beg; i<number_of_particles; i+=number_of_cores){ //возможно тут стоит создать новую переменную, чтобы все подряд не обращались к number_of_particles
			if ((rude_coords[i][0] >= 1) and (rude_coords[i][1] >= 1) and (rude_coords[i][2] >= 1) 
			and (rude_coords[i][0] < field_size-1) and (rude_coords[i][1] < field_size-1) and (rude_coords[i][2] < field_size-1)){
				//if (field[rude_coords[i][0]][rude_coords[i][1]][rude_coords[i][2]].front() != i){
					//std::cout << "SOS!" << field[rude_coords[i][0]][rude_coords[i][1]][rude_coords[i][2]].front() << '\n';
				//}
				for (long int x = -1; x<2; x++){
					for (long int y = -1; y<2; y++){
						for (long int z = -1; z<2; z++){
							std::list <long int> :: iterator cell;
							std::list <long int> :: iterator cell_end = field[rude_coords[i][0]+x][rude_coords[i][1]+y][rude_coords[i][2]+z].end();
							for (cell = field[rude_coords[i][0]+x][rude_coords[i][1]+y][rude_coords[i][2]+z].begin();
							    cell != cell_end; cell++){
								//print("F",i,j)
								if (i < *cell){
									double R6 = std::pow(sum(power(cur_coords[i]-cur_coords[*cell], 2)),3);
									local_energies += (1/(2*R6)-1)/R6 - U_piece;
								}
							}
						}
					}
				}
			}
		}

	}

	/*
	def check_energies(self):
		#нахождение точек со скачками по энергии
		for i in range (0, max_step-1):
			if (energies[i+1]-energies[i]) >1:
				print(i, '\t', energies[i], ' ', energies[i+1], '\n', data[i], '\n')

	*/
	void loading(){
		double pr = step*1000/max_step/10.0;
		if (percent < pr){
			percent = pr;
			//std::cout << pr << "%\t" << std::time(0) - time << " E = " << energy(step) << "s\n";
			std::cout << pr << "%\t" << std::time(0) - time << "s\n";
		}
	}
};

/*
class Visual:

	def __init__(self):
		body = Calc(125, False)
		body_calc_status = False


	def show_energy(self):
		if not body_calc_status:
			body.evolve()
			body_calc_status = True

		energies = np.zeros(shape = (body.max_step-1))
		for i in range (0, body.max_step-1):
			energies[i] = body.energy(i+1)[2] 

		plt.plot(np.arange(body.max_step-1), energies)
		plt.show()



	def status_body(self):

		if not body_calc_status:
			body.evolve()
			body_calc_status = True

		#T = time.clock()
		size = 3*100*body.R_intr
		ranges = np.zeros(shape = (size))
		quantity_of_averaging = min(body.max_step, 1001)
		center_num = int(body.number_of_particles/2 + body.number_of_particles**(2/3)/2 + body.number_of_particles**(1/3)/2)
		print(body.data[body.max_step-1][center_num])
		for i in range (1, quantity_of_averaging):
			data = body.data[body.max_step-i]
			for j in range (0, body.number_of_particles):
				r = int(100*np.power(np.sum(np.power(data[j]- data[center_num], 2)), 0.5))
				if r < size:
					ranges[r] += 1
			#print(f'{100*i//quantity_of_averaging}%\t{round(time.clock()-T)}s')

		plt.plot(np.arange(size)/100, ranges)
		plt.show()		


	def output(self):
		if not body_calc_status:
			body.evolve()
			body_calc_status = True
		body.output()


	def wall_output(self):
		wall = Calc(1, True)
		wall.output_wall()
*/

/*
extern "C" {
    Calc* Calc_new(){ return new Calc(125); }
    void Calc_evolve(Calc* c){ c->slow_evolve(); }
}
*/

int main(){
	/*
	for (long int i = 1; i < 16; i++){
		Calc c = Calc(i*i*i);
		//c.output_wall();

		c.multicore_status = false;
		time_t time_evolve1 = std::time(0);
		c.evolve();
		time_evolve1 = std::time(0) - time_evolve1;

		c.multicore_status = true;
		time_t time_evolve2 = std::time(0);
		c.evolve();
		time_evolve2 = std::time(0) - time_evolve2;

		std::cout << "particles: " << i*i*i << " multicores = " << false << ' ' << time_evolve1 << "s \n";
		std::cout << "particles: " << i*i*i << " multicores = " << true << ' ' << time_evolve2 << "s \n";
	}
	*/
	
	Calc c = Calc(343);
	c.wall_status = true;
	c.saving_status = false;
	c.calculus("VelCM");
	//c.output_wall();
	//c.end_analysis();
}

