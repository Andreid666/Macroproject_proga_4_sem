#include <boost/compute.hpp>
#include <CL>
using namespace boost;

/*
std::vector<float> std_vector(10);
compute::vector<float> compute_vector(std_vector.begin(), std_vector.end(), queue); 



compute::copy(compute_vector.begin(), compute_vector.end(), std_vector.begin(), queue);


auto device = compute::system::default_device(); // устройство по умолчанию это видеокарта
auto context = compute::context::context(device); // обычное объявление переменной
auto queue = compute::command_queue(context, device); // аналогично к предыдущему
*/


int main(){
	auto device = compute::system::default_device(); // устройство по умолчанию это видеокарта
	std::cout << device.name() << std::endl; 
}