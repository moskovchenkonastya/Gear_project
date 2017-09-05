#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <vector>
#include <tuple>
#include <memory>
#include <math.h>
#include <string>

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::vector;
using std::tuple;
using std::tie;
using std::make_tuple;
using std::shared_ptr;

#include "io.h"
#include "matrix.h"
#include "MyObject.h"

#define WHITE_COLOR make_tuple(255,255,255)
#define BLACK_COLOR make_tuple(0,0,0)
#define DYN_COLOR(color) make_tuple(0,color,0)

#define MIN(x,y) (x<y?x:y)
#define MAX(x,y) (x<y?y:x)

#define INF 1E20

typedef Matrix<uint> BinMatrix;

typedef Matrix<float> DtMatrix;


uint get_lim(const Image& in) {

    int max_h = 0;
   
    uint i, j;
     uint r, g, b;


    for (i =0 ; i< in.n_rows; i++) {
        for(j =0 ; j< in.n_cols; j++) {
            tie(r, g, b) = in(i, j);
            int bright = 0.3 * r + 0.59* g + 0.11 * b; // узнали яркость по RGB    
            max_h = (max_h> bright ? max_h : bright);
        }  
    }

    return max_h;


}


BinMatrix binarization(const Image& in) 
{
    BinMatrix bin (in.n_rows, in.n_cols);

    uint i, j;

    uint r, g, b, lim = get_lim(in);


    for (i =0 ; i< in.n_rows; i++) {
        for(j =0 ; j< in.n_cols; j++) {
            tie(r, g, b) = in(i, j);
            int bright = 0.3 * r + 0.59 * g + 0.11 * b; // узнали яркость по RGB

            bin(i,j) = ( bright > lim * 0.35 ? 1 : 0);
 
        }

    }

    return bin;
}


Image convert_matrix_to_image (const BinMatrix& in) 
{
    
    Image res (in.n_rows, in.n_cols);

    for (uint i = 0; i < in.n_rows; ++i) {
        for (uint j = 0; j < in.n_cols; ++j) {
            if (in(i,j) == 0 )
                res(i,j) = BLACK_COLOR;
            else
                res(i,j) = DYN_COLOR (in(i,j) * 40);
        }
    }

    return res;

}

// взято из лекции
void Fill(const BinMatrix& in, BinMatrix& labels, uint i, uint j, int L) 
{

    if( (labels(i, j) == 0) & (in(i,j) == 1)) 
    {
        labels(i, j) = L;

        if( i > 1 )
            Fill(in, labels, i - 1, j, L);

        if (i < in.n_rows - 1)
            Fill(in, labels, i + 1, j, L);
        
        if (j > 1)
           Fill(in, labels, i, j - 1, L);
        
        if (j < in.n_cols - 1)
           Fill(in, labels, i, j + 1, L);
    }
}

void segmentation(const BinMatrix& in, BinMatrix& labels) 
{
    uint i, j, L = 1;

    for (i = 0 ; i < in.n_rows; i++) 
        for(j = 0 ; j < in.n_cols; j++) {
            if (in(i,j) == 1 && labels(i, j) == 0) {
                // функция из лекции
                Fill(in, labels, i, j, L++);


            }
       }
}


tuple<uint,uint,uint,uint> 
get_segment(const BinMatrix& in, const uint L)
{
    uint x1 = in.n_rows, x2 = 0;
    uint y1 = in.n_cols, y2 = 0;
    uint i,j;
    for (i =0 ; i < in.n_rows; i++)
        for(j =0 ; j < in.n_cols; j++){

            if (in(i,j) == L) {

                x1 = MIN (x1, i);
                x2 = MAX (x2, i);

                y1 = MIN (y1, j);
                y2 = MAX (y2, j);

            }    
    }

 /*
    printf("%d:+ %4d %4d %4d %4d \n", L, x1, y1, x2,y2);
    printf("%d:- %4d %4d %4d %4d \n", L, x1, y1, x2-x1,y2-y1);
*/

    return make_tuple (x1, y1, x2-x1,y2-y1);

}

uint get_segment_counts(const BinMatrix& in)
{
    uint i, j, segments = 0;
     for (i = 0 ; i < in.n_rows; i++)
        for(j = 0 ; j < in.n_cols; j++)    
            segments = MAX (segments, in(i,j));

    return segments;
}

tuple<uint, uint, uint, uint>
get_radius (const BinMatrix& in, const DtMatrix& dt_map, const uint L) 
{
    
    uint i,j,x, cx, cy;
    
    uint min_r = 0;

    
    for(i = 0; i < in.n_rows; i++){
        for(j = 0; j < in.n_cols; j++) {

            if(in(i, j) == L){
                x = dt_map(i,j);
                if(min_r < x){
                    min_r = x;
                    cx = i;
                    cy = j;
                }
                
            }
        }
    }

    uint max_dist = 0, dist;

    for (i = 1 ; i < in.n_rows ; i++){
        for(j = 1 ; j < in.n_cols; j++) {

            if (in(i,j) == L) {
                dist = sqrt ((i - cx) * (i - cx) + (j-cy) * (j-cy) );
                max_dist = MAX(max_dist, dist);
            }
        }
    }
    
    return make_tuple(cy, cx, min_r, max_dist);
}



float *dt1(float *f, int n) {

    //cout << "dt1"<< endl;
    float *d = new float[n];
    int *v = new int[n];
    float *z = new float[n + 1];
    int k = 0;
    v[0] = 0;
    z[0] = -INF;
    z[1] = +INF;
    for (int q = 1; q <= n - 1; q++) {
            float s  = ((f[q]+q * q)-(f[v[k]]+v[k] * v[k])) / (2 * q-2*v[k]);
            while (s <= z[k]) {
                k--;
                s  = ((f[q]+ q * q)-(f[v[k]]+(v[k] * v[k])))/(2*q-2*v[k]);
        }
        k++;
        v[k] = q;
        z[k] = s;
        z[k + 1] = +INF;
  }

  k = 0;
  for (int q = 0; q <= n -1; q++) {
        while (z[k+1] < q)
        k++;
        d[q] = (q-v[k])* (q-v[k]) + f[v[k]];
  }

  delete [] v;
  delete [] z;
  return d;
}

// неправильно работает определение центра!

void dt2(DtMatrix in, uint h, uint w) 
{
    uint i, j;
    float *f = new float[MAX(h,w)];
    
  
    for (i = 0; i < w; i++) {
        for (j = 0; j < h; j++) {
            f[j] = in(j,i);
        }


        float *d = dt1(f, h);
        for (j = 0; j < h; j++){ 
            in(j,i) = d[j];
        }
        delete [] d;
    }

   
    for (j = 0; j < h; j++) {
        for (i = 0; i < w; i++) {
            f[i] = in(j,i);
        }   
        float *d = dt1(f, w);
        for (i = 0; i < w; i++) {
            in(j,i) = d[i];
        }
        delete [] d;
    }
    delete f;
}

DtMatrix dt3(BinMatrix in, uint h, uint w) {

    
    DtMatrix out(h, w);
    for (uint i = 0; i < h; i++) {
        for (uint j = 0; j < w; j++) {
            if (in(i,j) == 1)
                out(i,j) = INF;
            else
                out(i,j) = 0;
    }
  }
  dt2(out, h, w);
  return out;
}


shared_ptr<IObject> get_gear_info(BinMatrix in, DtMatrix dt_map, const uint L)
{
    uint left,top,w,h;
    tie (left,top,w,h)  = get_segment(in,L);

    // посчитаем площадь

    uint i, j, x_cm = 0, y_cm = 0, N = 0, cx, cy;
    for (i = 0; i < in.n_rows; i++){ 
        for(j = 0; j < in.n_cols; j++) {
            if (in(i, j) == L){
                N++;
                x_cm += i;
                y_cm += j;
            }
        }
    }

    BinMatrix dt(in.n_rows, in.n_cols);
    
    uint in_rad, out_rad;
    tie(cx, cy,in_rad,out_rad) = get_radius (in, dt_map, L);
    tuple<int, int> in_location = make_tuple(cx, cy);
    // является ли сломаной

    // количество зубцов

    shared_ptr<IObject> ret (new Gear(in_location, in_rad, out_rad));

    return  ret;
}


shared_ptr<IObject> get_axis_info(const BinMatrix& in, const uint L)
{
    uint left,top,w,h;
    tie (left,top,w,h)  = get_segment(in,L);

    tuple<int, int> in_location = make_tuple (top + h/2, left + w/2);

    shared_ptr<IObject> ret (new Axis(in_location));
    return  ret;
}

int is_axis(const BinMatrix& in, const DtMatrix& dt_map, const uint L) {

    uint left,top,w,h;
    tie (left,top,w,h)  = get_segment(in,L);

    uint in_rad, out_rad, cx,cy;
    tie(cx, cy,in_rad, out_rad) = get_radius (in, dt_map, L);

    if (in_rad - out_rad < 0.005)
        return true;

    return false;
}



tuple<int, vector<shared_ptr<IObject>>, Image>
repair_mechanism(const Image& in, const char* argv)
{
    // Base: return array of found objects and index of the correct gear
    // Bonus: return addit ional parameters of gears
    auto object_array = vector<shared_ptr<IObject>>();
    int result_idx = 0;

    // бинаризация
    BinMatrix binary = binarization(in);

    BinMatrix labels (in.n_rows, in.n_cols); 
     
    segmentation(binary, labels);


    DtMatrix dt_map(binary.n_rows, binary.n_cols);
    

    dt_map = dt3(binary, binary.n_rows, binary.n_cols);
    Image dst_image(binary.n_rows,binary.n_cols);
    
    for (uint i = 0; i < binary.n_rows; i++) {
        for (uint j = 0; j < binary.n_cols; j++) {
            dt_map(i,j) = sqrt(dt_map(i,j));
            //dst_image(i, j) = make_tuple(dt_map(i,j), dt_map(i,j), dt_map(i,j));
        }    
    }

    
    
    for (uint i = 1; i <= get_segment_counts(labels); i++) {

        if (is_axis(labels, dt_map, i)) {
            object_array.push_back(get_axis_info (labels, i));
        } else{
            object_array.push_back(get_gear_info (labels, dt_map, i));
        }       
    }

    // берем адрес картинки, подается на вход
    // от нее вычитаем 4 символа
    
    string answer = argv.c_str - 4;

    string variant[3] = {"_1.bmp", "_2.bmp", "_3.bmp"};


    for (uint i = 0; i < 3; i++){

        Image var_image = load_image(argv + variant[i]);
        binary = binarization(var_image);
        



    }
    


    //for(uint i = 0; i < in )
    // сохраняем полученную матрицу в файл
    //BinMatrix bord (in.n_rows, in.n_cols); 
    dst_image = convert_matrix_to_image (labels);
   // Image dst_image = convert_matrix_to_image (binary);

    return make_tuple(result_idx, object_array, dst_image);
}



int main(int argc, char **argv)
{
    if (argc != 4)
    {
        cout << "Usage: " << endl << argv[0]
             << " <in_image.bmp> <out_image.bmp> <out_result.txt>" << endl;
        return 0;
    }

    try {
        Image src_image = load_image(argv[1]);
        ofstream fout(argv[3]);

        vector<shared_ptr<IObject>> object_array;
        Image dst_image;
        int result_idx;
        tie(result_idx, object_array, dst_image) = repair_mechanism(src_image, argv[1]);
        save_image(dst_image, argv[2]);

        fout << result_idx << endl;
        fout << object_array.size() << endl;
        for (const auto &obj : object_array)
            obj->Write(fout);

    } catch (const string &s) {
        cerr << "Error: " << s << endl;
        return 1;
    }
}
