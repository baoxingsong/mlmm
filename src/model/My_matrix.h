//
// Created by baoxing on 2/21/18.
//

#ifndef MLMM_CPP_MY_MATRIX_H
#define MLMM_CPP_MY_MATRIX_H
#include <sstream>
#include <cstring>
template <typename T>

class My_matrix {
    private:
        T ** matrix;
        uint64_t num_row;
        uint64_t num_column;
    public:
        My_matrix( const uint64_t & _num_row, const uint64_t & _num_column );
        My_matrix( const My_matrix<T> & obj);
        My_matrix( const My_matrix<T> && obj);
        void set_values_zero();
        My_matrix();
        ~My_matrix();
        T ** get_matrix() const;
        const uint64_t & get_num_row() const;
        const uint64_t & get_num_column() const;
        My_matrix &operator=(const My_matrix &);
        void value_copy(const My_matrix<T> & obj);
};


template <typename T>
My_matrix<T>::My_matrix(){
    this->num_row=0;
    this->num_column=0;
    this->matrix = new T* [0];
}

template <typename T>
My_matrix<T>::My_matrix( const uint64_t & _num_row, const uint64_t & _num_column ){
    this->num_row=_num_row;
    this->num_column=_num_column;
    this->matrix = new T* [num_row];
    for( int i=0; i<num_row; ++i ){
        this->matrix[i] = new T[num_column];
    }
}

template <typename T>
My_matrix<T>::My_matrix( const My_matrix<T> & obj){
    this->num_row=obj.get_num_row();
    this->num_column=obj.get_num_column();
    this->matrix = new T* [num_row];
    for( int i=0; i<num_row; ++i ){
        this->matrix[i] = new T[num_column];
        memcpy(matrix[i], obj.get_matrix()[i], sizeof(T)*num_column );
    }
}

template <typename T>
void My_matrix<T>::value_copy(const My_matrix<T> & obj){
    for( int i=0; i<this->num_row; ++i ){
        memcpy(matrix[i], obj.get_matrix()[i], sizeof(T)*this->num_column );
    }
}
template <typename T>
My_matrix<T>::My_matrix( const My_matrix<T> && obj){
    this->num_row=obj.get_num_row();
    this->num_column=obj.get_num_column();
    this->matrix = new T* [num_row];
    for( int i=0; i<num_row; ++i ){
        this->matrix[i] = new T[num_column];
        memcpy(matrix[i], obj.get_matrix()[i], sizeof(T)*num_column );
    }
}

template <typename T>
void  My_matrix<T>::set_values_zero(){
    int i;
    int j;
    T value=0;
    for( i=0; i<this->num_row; ++i ){
        for( j=0; j<this->num_column; ++j ){
            this->matrix[i][j] = value;
        }
    }

}

template <typename T>
My_matrix<T> & My_matrix<T>::operator=(const My_matrix<T> & obj){
    this->num_row=obj.get_num_row();
    this->num_column=obj.get_num_column();
    this->matrix = new T* [this->num_row];
    for( int i=0; i<this->num_row; ++i ){
        this->matrix[i] = new T[this->num_column];
        memcpy(matrix[i], obj.get_matrix()[i], sizeof(T)*this->num_column );
    }
    return * this;
}


template <typename T>
My_matrix<T>::~My_matrix(){
    for ( int i=0; i<num_row; ++i ){
        delete [] this->matrix[i];
    }
    delete [] this->matrix;
}



template <typename T>
T ** My_matrix<T>::get_matrix() const{
    return this->matrix;
}

template <typename T>
const uint64_t & My_matrix<T>::get_num_row() const{
    return this->num_row;
}

template <typename T>
const uint64_t & My_matrix<T>::get_num_column() const{
    return this->num_column;
}

#endif //MLMM_CPP_MY_MATRIX_H
