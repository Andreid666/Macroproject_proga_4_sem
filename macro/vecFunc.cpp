#include <vector>
#include <iostream>
#include <ctime>
#include <random>
#include <cmath>

template <typename T>
std::vector<T> operator + (const std::vector<T> & v1, const std::vector<T> & v2){
	long int size = v1.size();
	if (size == v2.size()){
		std::vector<T> sum;
		sum.resize(size);
		for (long int i =0; i < size; i++){
			sum[i] = v1[i] + v2[i];
		}
		return sum;
	}
	else {
		std::cout << "xxx";
	}
}

template <typename T>
std::vector<T> operator - (const std::vector<T> & v1, const std::vector<T> & v2){
	long int size = v1.size();
	if (size == v2.size()){
		std::vector<T> sum;
		sum.resize(size);
		for (long int i =0; i < size; i++){
			sum[i] = v1[i] - v2[i];
		}
		return sum;
	}
	else {
		std::cout << "xxx";
	}
}

template <typename T>
T operator * (const std::vector<T> & v1, const std::vector<T> & v2){
	long int size = v1.size();
	if (size == v2.size()){
		T sum = 0;
		for (long int i =0; i < size; i++){
			sum += v1[i]*v2[i];
		}
		return sum;
	}
	else {
		std::cout << "xxx";
	}
}

template <typename T1, typename T2>
std::vector<T1> operator * (const std::vector<T1> & v, const T2 c){
	long int size = v.size();
	std::vector<T1> prod;
	prod.resize(size);
	for (long int i =0; i < size; i++){
		prod[i] = c*v[i];
	}
	return prod;
}

template <typename T1, typename T2>
std::vector<T1> operator * (const T2 c, const std::vector<T1> & v){
	long int size = v.size();
	std::vector<T1> prod;
	prod.resize(size);
	for (long int i =0; i < size; i++){
		prod[i] = c*v[i];
	}
	return prod;
}

template <typename T1, typename T2>
std::vector<T1> operator + (const std::vector<T1> & v, const T2 c){
	long int size = v.size();
	std::vector<T1> prod;
	prod.resize(size);
	for (long int i =0; i < size; i++){
		prod[i] = c+v[i];
	}
	return prod;
}

template <typename T1, typename T2>
std::vector<T1> operator + (const T2 c, const std::vector<T1> & v){
	long int size = v.size();
	std::vector<T1> prod;
	prod.resize(size);
	for (long int i =0; i < size; i++){
		prod[i] = c+v[i];
	}
	return prod;
}

template <typename T1, typename T2>
std::vector<T1> operator - (const std::vector<T1> & v, const T2 c){
	long int size = v.size();
	std::vector<T1> prod;
	prod.resize(size);
	for (long int i =0; i < size; i++){
		prod[i] = v[i]-c;
	}
	return prod;
}

template <typename T1, typename T2>
std::vector<T1> operator / (const std::vector<T1> & v, const T2 c){
	long int size = v.size();
	std::vector<T1> del;
	del.resize(size);
	for (long int i =0; i < size; i++){
		del[i] = (T1) (v[i]/c);
	}
	return del;
}

template <typename T>
std::vector<T> operator += (std::vector<T> & v1, const std::vector<T> & v2){
	long int size = v1.size();
	if (size == v2.size()){
		for (long int i =0; i < size; i++){
			v1[i] += v2[i];
		}
		return v1;
	}
	else {
		std::cout << "xxx";
	}
}

template <typename T1, typename T2>
std::vector<T1> operator += (std::vector<T1> & v, const T2 c){
	long int size = v.size();
		for (long int i =0; i < size; i++){
			v[i] += c;
		}
		return v;
}

template <typename T>
std::vector<T> operator -= (std::vector<T> & v1, const std::vector<T> & v2){
	long int size = v1.size();
	if (size == v2.size()){
		for (long int i =0; i < size; i++){
			v1[i] -= v2[i];
		}
		return v1;
	}
	else {
		std::cout << "xxx";
	}
}

template <typename T1, typename T2>
std::vector<T1> operator -= (std::vector<T1> & v, const T2 c){
	long int size = v.size();
		for (long int i =0; i < size; i++){
			v[i] -= c;
		}
		return v;
}

//Заполнение вектора рандомными числами от 0 до 1
void random_v(std::vector<double> & v){
	long int size = v.size();
	std::mt19937 gen(clock());
	double transl = -gen.min();
	double norm = 1.0/(gen.max() - gen.min());
	for (long int i=0; i<size; i++){
		v[i] = norm*gen() + norm;
	}
}

template <typename T>
std::vector<std::vector<T>> operator + (const std::vector<std::vector<T>> & v1, const std::vector<std::vector<T>> & v2){
	long int size = v1.size();
	if (size == v2.size()){
		std::vector<std::vector<T>> sum;
		sum.resize(size);
		for (long int i =0; i < size; i++){
			sum[i] = v1[i] + v2[i];
		}
		return sum;
	}
	else {
		std::cout << "xxx";
	}
}

template <typename T>
std::vector<std::vector<T>> operator - (const std::vector<std::vector<T>> & v1, const std::vector<std::vector<T>> & v2){
	long int size = v1.size();
	if (size == v2.size()){
		std::vector<std::vector<T>> sum;
		sum.resize(size);
		for (long int i =0; i < size; i++){
			sum[i] = v1[i] - v2[i];
		}
		return sum;
	}
	else {
		std::cout << "xxx";
	}
}

template <typename T>
T operator * (const std::vector<std::vector<T>> & v1, const std::vector<std::vector<T>> & v2){
	long int size = v1.size();
	if (size == v2.size()){
		T sum = 0;
		for (long int i =0; i < size; i++){
			sum += v1[i]*v2[i];
		}
		return sum;
	}
	else {
		std::cout << "xxx";
	}
}

template <typename T1, typename T2>
std::vector<std::vector<T1>> operator * (const std::vector<std::vector<T1>> & v, const T2 c){
	long int size = v.size();
	std::vector<std::vector<T1>> prod;
	prod.resize(size);
	for (long int i =0; i < size; i++){
		prod[i] = c*v[i];
	}
	return prod;
}

template <typename T1, typename T2>
std::vector<std::vector<T1>> operator * (const T2 c, const std::vector<std::vector<T1>> & v){
	long int size = v.size();
	std::vector<std::vector<T1>> prod;
	prod.resize(size);
	for (long int i =0; i < size; i++){
		prod[i] = c*v[i];
	}
	return prod;
}

template <typename T1, typename T2>
std::vector<std::vector<T1>> operator + (const std::vector<std::vector<T1>> & v, const T2 c){
	long int size = v.size();
	std::vector<std::vector<T1>> prod;
	prod.resize(size);
	for (long int i =0; i < size; i++){
		prod[i] = c+v[i];
	}
	return prod;
}

template <typename T1, typename T2>
std::vector<std::vector<T1>> operator + (const T2 c, const std::vector<std::vector<T1>> & v){
	long int size = v.size();
	std::vector<std::vector<T1>> prod;
	prod.resize(size);
	for (long int i =0; i < size; i++){
		prod[i] = c+v[i];
	}
	return prod;
}


template <typename T1, typename T2>
std::vector<std::vector<T1>> operator - (const std::vector<std::vector<T1>> & v, const T2 c){
	long int size = v.size();
	std::vector<std::vector<T1>> prod;
	prod.resize(size);
	for (long int i =0; i < size; i++){
		prod[i] = v[i]-c;
	}
	return prod;
}

template <typename T1, typename T2>
std::vector<std::vector<T1>> operator / (const std::vector<std::vector<T1>> & v, const T2 c){
	long int size = v.size();
	std::vector<std::vector<T1>> del;
	del.resize(size);
	for (long int i =0; i < size; i++){
		del[i] = v[i]/c;
	}
	return del;
}

template <typename T>
std::vector<std::vector<T>> operator += (std::vector<std::vector<T>> & v1, const std::vector<std::vector<T>> & v2){
	long int size = v1.size();
	if (size == v2.size()){
		for (long int i =0; i < size; i++){
			v1[i] += v2[i];
		}
		return v1;
	}
	else {
		std::cout << "xxx";
	}
}

template <typename T1, typename T2>
std::vector<std::vector<T1>> operator += (std::vector<std::vector<T1>> & v, const T2 c){
	long int size = v.size();
		for (long int i =0; i < size; i++){
			v[i] += c;
		}
		return v;
}

template <typename T>
std::vector<std::vector<T>> operator -= (std::vector<std::vector<T>> & v1, const std::vector<std::vector<T>> & v2){
	long int size = v1.size();
	if (size == v2.size()){
		for (long int i =0; i < size; i++){
			v1[i] -= v2[i];
		}
		return v1;
	}
	else {
		std::cout << "xxx";
	}
}

template <typename T1, typename T2>
std::vector<std::vector<T1>> operator -= (std::vector<std::vector<T1>> & v, const T2 c){
	long int size = v.size();
	for (long int i =0; i < size; i++){
		v[i] -= c;
	}
	return v;
}

template <typename T>
std::vector<long int> rounder(const std::vector<T> & v){
	long int size = v.size();
	std::vector<long int> res;
	res.resize(size);
	for (long int i =0; i < size; i++){
		res[i] = (long int) v[i];
	}
	return res;
}

template <typename T>
std::vector<std::vector<long int>> rounder(const std::vector<std::vector<T>> & v){
	long int size = v.size();
	std::vector<std::vector<long int>> res;
	res.resize(size);
	for (long int i =0; i < size; i++){
		res[i] = rounder(v[i]);
	}
	return res;
}

template <typename T>
T sum(const std::vector<T> & v){
	long int size = v.size();
	T summ = 0;
	for (int i=0; i < size; i++){
		summ += v[i];
	}
	return summ;
}

template <typename T>
T sum(const std::vector<std::vector<T>> & v){
	long int size = v.size();
	T summ = 0;
	for (int i=0; i < size; i++){
		summ += sum(v[i]);
	}
	return summ;
}

template <typename T1, typename T2>
std::vector<T1> power(const std::vector<T1> & v, T2 p){
	long int size = v.size();
	std::vector<T1> res;
	res.resize(size);
	for (int i=0; i < size; i++){
		res[i] = pow(v[i], p);
	}
	return res;
}

template <typename T1, typename T2>
std::vector<std::vector<T1>> power(const std::vector<std::vector<T1>> & v, T2 p){
	long int size = v.size();
	std::vector<std::vector<T1>> res;
	res.resize(size);
	for (int i=0; i < size; i++){
		res[i] = power(v[i], p);
	}
	return res;
}

