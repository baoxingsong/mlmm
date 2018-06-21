//
// Created by baoxing on 2/21/18.
//

#ifndef MLMM_CPP_MY_VECTOR_H
#define MLMM_CPP_MY_VECTOR_H

#include <sstream>
#include <vector>
#include <cstring>

template <typename T>
class My_Vector {
    private:
        T * array;
        size_t length=0;
    public:
        My_Vector(const size_t & _length);
        My_Vector(const size_t & _length, T * _array);
        My_Vector( const My_Vector & a);
        My_Vector( const My_Vector && a);
        My_Vector( const std::vector<T> & a);
        My_Vector();
        ~My_Vector();
        T * get_array() const;
//        void set_array(double * _array);
        const size_t & get_length() const;
        void reSet(const size_t & _length);
//        void set_length( const unsigned long & new_length );
        My_Vector<T> &operator=(const My_Vector<T> &);
        void set_values_zero();

};

template <typename T>
My_Vector<T>::My_Vector(const size_t & _length){
        this->length=_length;
        this->array = new T[_length];
}

template <typename T>
My_Vector<T>::My_Vector(const size_t & _length, T * _array){
        this->length=_length;
        this->array=_array;
}

template <typename T>
My_Vector<T>::My_Vector(const std::vector<T> & a){
        this->length = a.size();
        this->array = new T[this->length];
        for( int i=0; i<this->length; ++i ){
                this->array[i]=a[i];
        }
}

template <typename T>
My_Vector<T>::My_Vector(){
}

template <typename T>
My_Vector<T>::My_Vector( const My_Vector & a){
        this->length=a.get_length();
        this->array = new T[this->length];
        memcpy(this->array, a.get_array(), sizeof(T)*this->length );
}

template <typename T>
My_Vector<T>::My_Vector( const My_Vector && a){
        this->length=a.get_length();
        this->array = new T[this->length];
        memcpy(this->array, a.get_array(), sizeof(T)*this->length );
}

template <typename T>
My_Vector<T>::~My_Vector(){
        if( this->length > 0 ){
                delete [] array;
        }
}

template <typename T>
T * My_Vector<T>::get_array() const{
    return this->array;
}
//
//void My_Vector::set_array(double * _array){
//    this->array=_array;
//}
template <typename T>
const size_t & My_Vector<T>::get_length() const{
    return this->length;
}
template <typename T>
void My_Vector<T>::reSet(const size_t & _length){
    if( this->length > 0 ){
        delete [] array;
    }
    this->length=_length;
    this->array = new T[_length];
}
//void My_Vector::set_length( const unsigned long & new_length ){
//    this->length=new_length;
//}

template <typename T>
My_Vector<T> &My_Vector<T>::operator=(const My_Vector<T> & a){
        this->length=a.get_length();
        this->array = new T[this->length];
        memcpy(this->array, a.get_array(), sizeof(T)*this->length );
        return * this;
}

template <typename T>
void  My_Vector<T>::set_values_zero(){
        size_t i;
        for( i=0; i<this->get_length(); ++i ){
                this->get_array()[i]=0;
        }
}

#endif //MLMM_CPP_MY_VECTOR_H
