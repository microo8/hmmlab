%module(docstring = "Module for creating and manipulating HMMs") libhmm
%include "std_vector.i"
%include "std_map.i"
%include "std_string.i"
%allowexception;

%{
#define SWIG_FILE_WITH_INIT
#include "libhmm.h"
using namespace std;
%}

%feature("ref")   RCObj "$this->inc_ref_num();"
%feature("unref") RCObj "$this->dec_ref_num();"

%include "data_structures.h"

%define vec_list_wrap(postfix, T)
    %template(vector ## postfix) std::vector<T >;
    %template(List ## postfix) List<T >;
%enddef

vec_list_wrap(Int, int)
vec_list_wrap(Double, double)
vec_list_wrap(Model, Model*)
vec_list_wrap(State, State*)
vec_list_wrap(Stream, Stream*)
vec_list_wrap(Gaussian, Gaussian*)
vec_list_wrap(ListGaussian, List<Gaussian* >)
vec_list_wrap(Vector, Vector*)
vec_list_wrap(SVector, SVector*)
vec_list_wrap(ListDouble, List<double>*)
vec_list_wrap(ListVector, List<Vector*>)
vec_list_wrap(ListListDouble, List<List<double>* >)
vec_list_wrap(ListListVector, List<List<Vector*>* >)
vec_list_wrap(HMMLab_Object, HMMLab_Object*)

%define map_dict_wrap(postfix, T, U)
    %template(map ## postfix) std::map<T , U >;
    %template(Dict ## postfix) Dict<T , U >;
%enddef

map_dict_wrap(StringHMMLab_Object, std::string, HMMLab_Object * )

%include "data_structures.h"

%rename(assign) *::operator=;
%rename(__getitem__) *::operator[];
%rename(inc) *::operator++;

%init {
    init();
}

%include "data_structures.h"
%include "vmlib.h"
%include "libhmm.h"
