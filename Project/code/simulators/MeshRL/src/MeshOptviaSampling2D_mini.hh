#include <string>
#include <vector>
#include <cassert>
#include <iostream>
#include <fstream>
#include <time.h>
#include <sstream>

#define PI 3.14159265358979323846264338327950288
#define TwoPI 6.2831853071795862

//#define PeriodicPoints
//#define no_boundary_sifting
//#define disconnected_regions

#define debug


double cc;
int indx;
double Q[1220];

double**_samples (NULL), _tol, _angle_tol; //sample coordinates (x,y) + boundary_id + 1 for corner 0 otherwise, tolerance 
size_t _num_samples_input, _num_samples, ** _del(NULL), _max_del(19);//number of samples, connectivity
size_t*_list = NULL;//neightbouring samples sorted list

double**_circum_ex = NULL;//x,y,r of the neighbour not-to-be-distributed delaunay circumcircles exclusion regions
size_t _num_circum_ex;//number of neighbour not-to-be-distributed delaunay circumcircles 
double**_circum_in = NULL;//x,y,r of the neighbour not-to-be-formed delaunay circumcircles inclusion regions
size_t _num_circum_in;//number of neighbour not-to-be-formed delaunay circumcircles 
double**_line_in = NULL;
size_t _num_line_in;
size_t _max_active_cell, *_active_cells_i(NULL), *_active_cells_j(NULL), *_tmp_active_cells_i(NULL), *_tmp_active_cells_j(NULL);
double* _angles = NULL;//used for sorting
double _smoothness_factor, _smoothness_angle;
double**_xy_org = NULL;
size_t **_ghost_points = NULL;

double _min_angle_input, _max_angle_input;
size_t _num_obtuse_input(0);//total number of obtuse angles in the input 

size_t _numSifting(0), _numSiftingAtt(0), _numInject(0), _numInjectRep(0);
size_t _numSifting_failed(0), _numSiftingAtt_failed(0), _numInject_failed(0), _numInjectRep_failed(0);

size_t _current_num_non_obtuse;//current number of non-obtuse triangles 
bool _first_read = true;

double angle_three_points_stathead(double *xp1, double* xp2, double* xp3)
{
	double x1, x2, x3, y1, y2, y3;

	x1 = xp1[0]; y1 = xp1[1];
	x2 = xp2[0]; y2 = xp2[1];
	x3 = xp3[0]; y3 = xp3[1];

	double l1 = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
	double l2 = sqrt((x3 - x2)*(x3 - x2) + (y3 - y2)*(y3 - y2));


	double A = (x2 - x1)*(x2 - x3);
	double B = (y2 - y1)*(y2 - y3);
	double val = (A + B) / (l1*l2);
	if (val<-1.0){ val = -1.0; }
	if (val>1.0){ val = 1.0; }

	double Angle = acos(val);
	double pi = 3.1415926535897931;
	Angle = (Angle * 180) / pi;

	return Angle;
}
void get_max_min_angles(size_t _num_p, size_t**_delaunay, double**_xp, double&max, double&min, size_t FindApex(size_t ii, size_t jj, size_t skip1, size_t skip2, size_t skip3, size_t skip4))
{

	size_t ip, ip1, ip2, i;
	double angle;
	max = 0.0;
	min = 180.0;
	for (ip = 0; ip<_num_p; ip++){
		for (i = 1; i <= _delaunay[ip][0]; i++){
			ip1 = _delaunay[ip][i];
			ip2 = FindApex(ip, ip1, _num_p, _num_p, _num_p, _num_p);
			if (ip2 == _num_p){ continue; }

			angle = angle_three_points_stathead(_xp[ip1], _xp[ip], _xp[ip2]);
			if (angle>max){ max = angle; }
			if (angle<min){ min = angle; }
		}
	}
}
double DistPlotter(double x1, double y1, double z1, double x2, double y2, double z2)
{
	double dx, dy, dz;
	dx = x1 - x2;
	dy = y1 - y2;
	dz = z1 - z2;
	dx *= dx;
	dy *= dy;
	dz *= dz;

	return dx + dy + dz;

}
void mark_point(double xx, double yy, double** samples, double** dots, size_t num_points, size_t num_dots, double r)
{
	std::fstream file("disks.ps", std::ios::app);

	double scale_x, scale_y, scale;

	double xmin(10E-8), ymin(10E-8), xmax(-10E-8), ymax(-10E-8);

	for (size_t i = 0; i < num_points; i++)
	{
		if (samples[i][0] < xmin) xmin = samples[i][0];
		if (samples[i][1] < ymin) ymin = samples[i][1];
		if (samples[i][0] > xmax) xmax = samples[i][0];
		if (samples[i][1] > ymax) ymax = samples[i][1];
	}

	double Lx = (xmax - xmin);
	double Ly = (ymax - ymin);

	scale_x = 6.5 / Lx;
	scale_y = 9.0 / Ly;

	if (scale_x < scale_y)
	{
		scale = scale_x;
	}
	else
	{
		scale = scale_y;
	}




	file << (xx - xmin)*scale << " " << (yy - ymin)*scale << " " << r*scale << " dark_blue_dot" << std::endl;

	//double l = 0.15*r*cos(45.0*3.14 / 180.0);

	//file<< (xx - l -xmin)*scale<<" "<< (yy - l -ymin)*scale<<"	"<< (xx + l -xmin)*scale<<" "<< (yy + l -ymin)*scale<<"	r_seg"<<std::endl;
	//file<< (xx - l -xmin)*scale<<" "<< (yy + l -ymin)*scale<<"	"<< (xx + l -xmin)*scale<<" "<< (yy - l -ymin)*scale<<"	r_seg"<<std::endl;


}
inline size_t FindApex(size_t ip, size_t ip1, size_t skip1, size_t skip2, size_t skip3, size_t skip4, size_t**connectivity)
{
	size_t i, j;

	for (i = 1; i <= connectivity[ip][0]; i++){
		if (connectivity[ip][i] == skip1 || connectivity[ip][i] == skip2 || connectivity[ip][i] == skip3 || connectivity[ip][i] == skip4){ continue; }
		for (j = 1; j <= connectivity[ip1][0]; j++){
			if (connectivity[ip1][j] == skip1 || connectivity[ip1][j] == skip2 || connectivity[ip][i] == skip3 || connectivity[ip][i] == skip4){ continue; }
			if (connectivity[ip][i] == connectivity[ip1][j]){

				return connectivity[ip][i];
			}
		}
	}

	return skip1;
}
void plot_n_layers(size_t n_layers, size_t ip_main, double** samples, size_t** connectivity, size_t num_points,
	double r, bool plot_samples_indices, bool plot_disks, bool show_connectivity, bool hide_ghost_points,
	bool *skip_list, bool use_skip_list, size_t num_dots, double**dots, bool voronoi, double**vor,
	bool use_red_green, bool*cond, bool colorize)
{
	std::fstream file("disks.ps", std::ios::out);
	

	file << "%!PS-Adobe-3.0" << std::endl;
	file << "72 72 scale     % one unit = one inch" << std::endl;


	file << "/quad_white      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " closepath" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " 1.0 1.0 1.0 setrgbcolor" << std::endl;
	file << "fill" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/orange_disk      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 1 0.72 0 setrgbcolor" << std::endl;
	//file << " fill" << std::endl;
	file << " 0.005 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/blue_disk      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0 0.33 1.0 setrgbcolor" << std::endl;
	//	file << " fill" << std::endl;
	file << " 0.005 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/orange_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 1 0.72 0 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/purple_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0.98 0 0.717 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/green_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0.25 0.75 0.25 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/quad_purp      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " closepath" << std::endl;
	file << "0.49 0.15 0.8 setrgbcolor" << std::endl;
	//	file << "fill" << std::endl;
	file << " 0.005 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/quad_bold      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " closepath" << std::endl;
	file << "0 0 0 setrgbcolor" << std::endl;
	//	file << "fill" << std::endl;
	file << " 0.01 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0 0 0 setrgbcolor" << std::endl;
	//file << " fill" << std::endl;
	file << " 0.006 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/black_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0 0 0 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.006 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/red_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 1 0 0 setrgbcolor" << std::endl;
	//	file <<47/255.0<<"	"<<79/255.0<<"	"<< 79/255.0<<" setrgbcolor" << std::endl;	
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/green_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << "0 0.392156862745098 0 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/orange_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 1 0.72 0 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/purple_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0.98 0 0.717 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/gray_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0.83 0.83 0.83 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/disk      % stack: x y r" << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0 0 0 setrgbcolor" << std::endl;
	//	file << " fill" << std::endl;
	file << " 0.003 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;


	file << "/r_seg      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " 0.005 setlinewidth" << std::endl;
	file << " 0 0 0 setrgbcolor" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/orange_seg      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " 0.005 setlinewidth" << std::endl;
	file << " 1 0.72 0 setrgbcolor" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/faint_seg      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " 0.0002 setlinewidth" << std::endl;
	file << " 1.0 0.67 0.67 setrgbcolor" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/g_seg      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " 0.005 setlinewidth" << std::endl;
	file << " 0.33 0.33 0.33 setrgbcolor" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;


	file << "/blue_tri      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " closepath" << std::endl;
	file << "0 0 1 setrgbcolor" << std::endl;
	//	file << "fill" << endl;
	file << " 0.002 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/red_tri      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " closepath" << std::endl;
	file << "1 0 0 setrgbcolor" << std::endl;
	file << "fill" << std::endl;
	file << " 0.002 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/green_tri      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " closepath" << std::endl;
	file << "0 0.392156862745098 0 setrgbcolor" << std::endl;
	//	file << "fill" << endl;
	file << " 0.002 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;


	double scale_x, scale_y, scale, shift_x, shift_y;

	double ss = r / sqrt(2.0);

	double xmin(samples[ip_main][0] - ss*n_layers), ymin(samples[ip_main][1] - ss*n_layers), xmax(samples[ip_main][0] + ss*n_layers), ymax(samples[ip_main][1] + ss*n_layers);

	double Lx = (xmax - xmin);
	double Ly = (ymax - ymin);

	scale_x = 6.5 / Lx;
	scale_y = 9.0 / Ly;

	if (scale_x < scale_y)
	{
		scale = scale_x;
		shift_x = 0.9;
		shift_y = (11.0 - (Lx*scale)) / 2.0;
	}
	else
	{
		scale = scale_y;
		shift_x = (8.3 - (Ly*scale)) / 2.0;
		shift_y = 1.35;
	}
	file << shift_x << " " << shift_y << " translate" << std::endl;
	if (num_dots > 0){
		for (size_t i = 0; i < num_dots; i++){
			if (dots[i][2]>1000){ continue; }
			file << (dots[i][0] - xmin)*scale << " " << (dots[i][1] - ymin)*scale << " " << sqrt(dots[i][2])*scale << " dot" << std::endl;
			file << (dots[i][0] - xmin)*scale << " " << (dots[i][1] - ymin)*scale << " " << sqrt(dots[i][2])*scale << " orange_disk" << std::endl;
			
		}
	}


	/*file << (0-xmin)* scale << "  " << (0-ymin)*scale << "  ";
	file << (0-xmin)*  scale << "  " << (1-ymin) * scale << "  ";
	file << (1-xmin)* scale << "  " << (1-ymin)* scale << "  ";
	file << (1-xmin)* scale << "  " << (0-ymin)*scale << "  ";
	file << "quad_bold"      << std::endl;*/




	/*if(plot_edges_list){
	for(size_t i=0;i<edges.size();i++){
	file<<(samples[edges[i][0]][0]-xmin)*scale<<" "<<(samples[edges[i][0]][1]-ymin)*scale<<" "<<(samples[edges[i][1]][0]-xmin)*scale<<" "<<(samples[edges[i][1]][1]-ymin)*scale<<" orange_seg"<<std::endl;
	}
	}*/

	if (plot_disks){
		for (size_t i = 0; i < num_points; i++){
			if (samples[i][0] < xmin || samples[i][1] < ymin || samples[i][0] > xmax || samples[i][1] > ymax){ continue; }
			if (use_skip_list && skip_list[i]) { continue; }

			//if (samples[i][2] < 0){ continue; }//removed

			//double r_i = SizingFunction(samples[i][0], samples[i][1], _sx, _sy, _ny, sizingFunc);
			//r_i = sqrt(r_i);

			file << (samples[i][0] - xmin)*scale << " " << (samples[i][1] - ymin)*scale << " " << r*scale / 30.2814 << " black_dot" << std::endl;
			file << (samples[i][0] - xmin)*scale << " " << (samples[i][1] - ymin)*scale << " " << r*scale << " disk" << std::endl;
			
		}
	}


	//colorize 
	if (colorize){
		for (size_t i = 0; i < num_points; i++){
			if (use_skip_list){
				if (skip_list[i]) continue;
			}
			size_t length = connectivity[i][0];
			//size_t j_plus = length;
			size_t ip_p = connectivity[i][length];

			for (size_t j = 1; j <= length; j++){

				size_t ip = connectivity[i][j];

				ip_p = FindApex(i, ip, num_points, num_points, num_points, num_points, connectivity);
				if (ip_p == num_points){
					continue;
				}

				if (angle_three_points_stathead(samples[ip], samples[i], samples[ip_p]) == 0){
					file << (samples[ip_p][0] - xmin)*scale << " " << (samples[ip_p][1] - ymin)*scale << " "
						<< (samples[ip][0] - xmin)*scale << " " << (samples[ip][1] - ymin)*scale << " "
						<< (samples[i][0] - xmin)*scale << " " << (samples[i][1] - ymin)*scale << " ";
					file << " red_tri" << std::endl;
				}
			}
		}
	}


	if (show_connectivity)
	{
		for (size_t i = 0; i < num_points; i++)
		{
			if (samples[i][0] < xmin || samples[i][1] < ymin || samples[i][0] > xmax || samples[i][1] > ymax) {
				continue;
			}

			//if (samples[i][0] < 0 || samples[i][1] < 0 || samples[i][0] > 1 || samples[i][1] > 1 )
			//	continue;

			if (use_skip_list && skip_list[i]) continue;


			int length = connectivity[i][0];
			int j;
			size_t ip_p = connectivity[i][length];
			size_t j_plus = length;
			for (j = 1; j <= length; j++)
			{
				size_t ip = connectivity[i][j];

				file << (samples[i][0] - xmin)*scale << " " << (samples[i][1] - ymin)*scale << " " << (samples[ip][0] - xmin)*scale << " " << (samples[ip][1] - ymin)*scale << " faint_seg" << std::endl;
				

				if (voronoi/*&&samples[i][2]==0*/){
					if (samples[i][2] == 0){
						file << (vor[i][j * 2 - 1] - xmin)*scale << " " << (vor[i][j * 2] - ymin)*scale << " " << (vor[i][j_plus * 2 - 1] - xmin)*scale << " " << (vor[i][j_plus * 2] - ymin)*scale << " g_seg" << std::endl;
						
						j_plus = j;
					}
					else{
						file << (vor[i][j * 2 - 1] - xmin)*scale << " " << (vor[i][j * 2] - ymin)*scale << " " << (vor[i][(j + 1) * 2 - 1] - xmin)*scale << " " << (vor[i][(j + 1) * 2] - ymin)*scale << " g_seg" << std::endl;
						
					}
				}
				continue;

				file << (samples[ip_p][0] - xmin)*scale << " " << (samples[ip_p][1] - ymin)*scale << " "
					<< (samples[ip][0] - xmin)*scale << " " << (samples[ip][1] - ymin)*scale << " "
					<< (samples[i][0] - xmin)*scale << " " << (samples[i][1] - ymin)*scale << " ";

				//		file<<" red_tri"<< std::endl;
				//		ip_p = ip;		
				//		continue;

				if (angle_three_points_stathead(samples[ip], samples[i], samples[ip_p]) == 0)
				{
					file << " red_tri" << std::endl;
				}
				if (angle_three_points_stathead(samples[ip], samples[i], samples[ip_p]) == 1)
				{
					file << " blue_tri" << std::endl;
				}
				if (angle_three_points_stathead(samples[ip], samples[i], samples[ip_p]) == 2)
				{
					file << " green_tri" << std::endl;
				}
				ip_p = ip;
				//				file<<(samples[i][0]-xmin)*scale<<" "<<(samples[i][1]-ymin)*scale<<" "<<(samples[ip][0]-xmin)*scale<<" "<<(samples[ip][1]-ymin)*scale<<" r_seg"<<std::endl;
			}
			if (voronoi /*&& samples[i][3]==0*/){
				if (samples[i][2] != 0){
					file << (vor[i][j * 2 - 1] - xmin)*scale << " " << (vor[i][j * 2] - ymin)*scale << " " << (samples[i][0] - xmin)*scale << " " << (samples[i][1] - ymin)*scale << " g_seg" << std::endl;
					file << (samples[i][0] - xmin)*scale << " " << (samples[i][1] - ymin)*scale << " " << (vor[i][(1) * 2 - 1] - xmin)*scale << " " << (vor[i][(1) * 2] - ymin)*scale << " g_seg" << std::endl;

				}
			}
		}
	}

	if (use_red_green){
		for (size_t i = 0; i < num_points; i++){
			if (samples[i][0] < xmin || samples[i][1] < ymin || samples[i][0] > xmax || samples[i][1] > ymax)
				continue;
			if (cond[i]){
				file << (samples[i][0] - xmin)*scale << " " << (samples[i][1] - ymin)*scale << " " << r*scale / 15.0 << " red_dot" << std::endl;
			}
			else{
				file << (samples[i][0] - xmin)*scale << " " << (samples[i][1] - ymin)*scale << " " << r*scale / 15.0 << " green_dot" << std::endl;
			}
		}
	}
	for (size_t i = 0; i<num_points; i++){
		if (samples[i][2] != 0){
			double min_ed_len(10E5), r_plot, dist;
			for (size_t V = 1; V <= connectivity[i][0]; V++){
				dist = DistPlotter(samples[i][0], samples[i][1], 0, samples[connectivity[i][V]][0], samples[connectivity[i][V]][1], 0);
				if (dist<min_ed_len){ min_ed_len = dist; }
			}
			r_plot = sqrt(min_ed_len) / 15;
			r_plot = 0.0005;
			if (samples[i][3] != 0){
				file << (samples[i][0] - xmin)*scale << " " << (samples[i][1] - ymin)*scale << " " << r_plot*scale << " purple_dot" << std::endl;
			}
			else if (samples[i][2] != 0){
				file << (samples[i][0] - xmin)*scale << " " << (samples[i][1] - ymin)*scale << " " << r_plot*scale << " orange_dot" << std::endl;

			}
		}
	}
	if (plot_samples_indices)
	{
		for (size_t i = 0; i < num_points; i++)
		{
			if (samples[i][0] < xmin || samples[i][1] < ymin || samples[i][0] > xmax || samples[i][1] > ymax){ continue; }

			if (use_skip_list && skip_list[i]){ continue; }

			if (samples[i][2] < 0){ continue; }//removed

			file << "/Times-Roman findfont" << std::endl;
			file << "0.06 scalefont" << std::endl;
			file << "0 0 0 setrgbcolor" << std::endl;
			file << "setfont" << std::endl;
			file << (samples[i][0] - xmin)*scale << " " << (samples[i][1] - ymin)*scale << " moveto" << std::endl;
			file << "(" << i;
			file << ") " << "show" << std::endl;
		}
	}

	if (hide_ghost_points)
	{
		file << (-10.0)* scale << "  " << (-10.0)*scale << "  ";
		file << (-10.0)*  scale << "  " << (0 - ymin) * scale << "  ";
		file << (10.0)* scale << "  " << (0 - ymin)* scale << "  ";
		file << (10.0)* scale << "  " << (-10)*scale << "  ";
		file << "quad_white" << std::endl;

		file << (-10.0)* scale << "  " << (-10.0)*scale << "  ";
		file << (-10.0)*  scale << "  " << (10.0) * scale << "  ";
		file << (-xmin)* scale << "  " << (10.0)* scale << "  ";
		file << (-xmin)* scale << "  " << (-10)*scale << "  ";
		file << "quad_white" << std::endl;

		file << (-10.0)* scale << "  " << (10.0)*scale << "  ";
		file << (10.0)*  scale << "  " << (10.0) * scale << "  ";
		file << (10.0)* scale << "  " << (1.0 - ymin)* scale << "  ";
		file << (-10.0)* scale << "  " << (1.0 - ymin)*scale << "  ";
		file << "quad_white" << std::endl;

		file << (1.0 - xmin)* scale << "  " << (-10.0)*scale << "  ";
		file << (1.0 - xmin)*  scale << "  " << (10.0) * scale << "  ";
		file << (10.0)* scale << "  " << (10)* scale << "  ";
		file << (10.0)* scale << "  " << (-10)*scale << "  ";
		file << "quad_white" << std::endl;
	}

}
void plot_dots_n_layers_orange_list(size_t n_layers, size_t ip_main, double** samples, size_t*list, double r)
{
	std::fstream file("disks.ps", std::ios::app);

	double scale_x, scale_y, scale;
	double ss = r / sqrt(2.0);

	double xmin(samples[ip_main][0] - ss*n_layers), ymin(samples[ip_main][1] - ss*n_layers), xmax(samples[ip_main][0] + ss*n_layers), ymax(samples[ip_main][1] + ss*n_layers);

	double Lx = (xmax - xmin);
	double Ly = (ymax - ymin);

	scale_x = 6.5 / Lx;
	scale_y = 9.0 / Ly;

	if (scale_x < scale_y){
		scale = scale_x;
	}
	else{
		scale = scale_y;
	}

	for (size_t i = 1; i <= list[0]; i++){
		file << (samples[list[i]][0] - xmin)*scale << " " << (samples[list[i]][1] - ymin)*scale << " " << r*scale / 20.2814 << " black_dot" << std::endl;
	}


}
void plot_single_dot_n_layers(size_t n_layers, size_t ip_main, double** samples, double xx, double yy, double r)
{
	std::fstream file("disks.ps", std::ios::app);

	double scale_x, scale_y, scale;

	double ss = r / sqrt(2.0);

	double xmin(samples[ip_main][0] - ss*n_layers), ymin(samples[ip_main][1] - ss*n_layers), xmax(samples[ip_main][0] + ss*n_layers), ymax(samples[ip_main][1] + ss*n_layers);

	double Lx = (xmax - xmin);
	double Ly = (ymax - ymin);

	scale_x = 6.5 / Lx;
	scale_y = 9.0 / Ly;

	if (scale_x < scale_y)
	{
		scale = scale_x;
	}
	else
	{
		scale = scale_y;
	}

	file << (xx - xmin)*scale << " " << (yy - ymin)*scale << " " << r*scale / 8 << " green_dot" << std::endl;




}
void plot_unit_box(const char*NodeFileName, double** samples, size_t** connectivity, size_t num_points, double r, bool colorize, bool plot_samples_indices, bool plot_disks, bool show_connectivity, bool hide_ghost_points, size_t skip1, size_t skip2, bool skip, bool *skip_list, bool use_skip_list, bool voronoi, double**vor)
{
	std::fstream file(NodeFileName, std::ios::out);
	
	
	file.precision(15);

	file << "%!PS-Adobe-3.0" << std::endl;
	file << "72 72 scale     % one unit = one inch" << std::endl;


	file << "/quad_white      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " closepath" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " 1.0 1.0 1.0 setrgbcolor" << std::endl;
	file << "fill" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/quad_bold      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " closepath" << std::endl;
	file << "0 0 0 setrgbcolor" << std::endl;
	//	file << "fill" << std::endl;
	file << " 0.01 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0 0 0 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.01 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/red_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 1 0 0 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.008 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/green_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << "0 0.392156862745098 0 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/orange_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 1 0.72 0 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/purple_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0.98 0 0.717 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/gray_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0.8 0.8 0.8 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.008 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/black_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0.0 0.0 0.0 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.008 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/dark_blue_dot      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0.0 0.0 1.0 setrgbcolor" << std::endl;
	file << " fill" << std::endl;
	file << " 0.008 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/disk      % stack: x y r" << std::endl;
	file << "{newpath" << std::endl;
	file << " 0 360 arc closepath" << std::endl;
	file << " 0 0 1 setrgbcolor" << std::endl;
	//	file << " fill" << std::endl;
	file << " 0.001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/r_seg      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " 1 0.8 0 setrgbcolor" << std::endl;
	file << " 0.00000001 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/orange_seg      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " 0.005 setlinewidth" << std::endl;
	file << " 1 0.72 0 setrgbcolor" << std::endl;
	file << " stroke" << std::endl;

	file << "} def" << std::endl;
	file << "/faint_seg      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " 0.0002 setlinewidth" << std::endl;
	file << " 1.0 0.67 0.67 setrgbcolor" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;


	file << "/g_seg      % stack: x  y " << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " 0.000001 setlinewidth" << std::endl;
	file << " 0.5 0.5 0.5 setrgbcolor" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;


	file << "/blue_tri      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " closepath" << std::endl;
	file << "0 0 1 setrgbcolor" << std::endl;
	//	file << "fill" << endl;
	file << " 0.002 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/red_tri      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " closepath" << std::endl;
	file << "1 0 0 setrgbcolor" << std::endl;
	file << "fill" << std::endl;
	file << " 0.0002 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/green_tri      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " closepath" << std::endl;
	file << "0 0.392156862745098 0 setrgbcolor" << std::endl;
	//	file << "fill" << endl;
	file << " 0.002 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;

	file << "/black_tri      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << std::endl;
	file << "{newpath" << std::endl;
	file << " moveto" << std::endl;
	file << " lineto" << std::endl;
	file << " lineto" << std::endl;
	file << " closepath" << std::endl;
	file << "0 0 0 setrgbcolor" << std::endl;
	//	file << "fill" << endl;
	file << " 0.0002 setlinewidth" << std::endl;
	file << " stroke" << std::endl;
	file << "} def" << std::endl;


	double scale_x, scale_y, scale, shift_x, shift_y;

	double xmin(10E-8), ymin(10E-8), xmax(-10E-8), ymax(-10E-8);

	for (size_t i = 0; i < num_points; i++)
	{
		if (samples[i][0] < xmin) xmin = samples[i][0];
		if (samples[i][1] < ymin) ymin = samples[i][1];
		if (samples[i][0] > xmax) xmax = samples[i][0];
		if (samples[i][1] > ymax) ymax = samples[i][1];
	}

	double Lx = (xmax - xmin);
	double Ly = (ymax - ymin);

	scale_x = 6.5 / Lx;
	scale_y = 9.0 / Ly;

	if (scale_x < scale_y)
	{
		scale = scale_x;
		shift_x = 0.9;
		shift_y = (11.0 - (Lx*scale)) / 2.0;
	}
	else
	{
		scale = scale_y;
		//shift_x = (8.3-(Ly*scale))/2.0;
		shift_x = 0.9;
		shift_y = 1.35;
	}

	file << shift_x << " " << shift_y << " translate" << std::endl;
	int j;


	/*
	file << (xmin + 3.0*sqrt(_r))* scale << "  " << (ymin + 3.0*sqrt(_r))*scale << "  ";
	file << (xmin + 3.0*sqrt(_r))*  scale << "  " << (ymax - 3.0*sqrt(_r))* scale << "  ";
	file << (xmax - 1.0*sqrt(_r))* scale << "  " << (ymax - 3.0*sqrt(_r))* scale << "  ";
	file << (xmax - 1.0*sqrt(_r))* scale << "  " << (ymin + 3.0*sqrt(_r))*scale << "  ";
	file << "quad_bold" << std::endl;*/

	/*file << (xmin)* scale << "  " << (ymin)*scale << "  ";
	file << (xmin)*  scale << "  " << (ymax) * scale << "  ";
	file << (xmax)* scale << "  " << (ymax)* scale << "  ";
	file << (xmax)* scale << "  " << (ymin)*scale << "  ";
	file << "quad_bold"      << std::endl;
	*/

	if (plot_disks){
		for (size_t i = 0; i < num_points; i++){
			if (use_skip_list){
				if (skip_list[i]) continue;
			}
			if ((i == skip1 || i == skip2) && skip) continue;
			file << (samples[i][0] - xmin)*scale << " " << (samples[i][1] - ymin)*scale << " " << r*scale / 30.2814 << " dot" << std::endl;
			file << (samples[i][0] - xmin)*scale << " " << (samples[i][1] - ymin)*scale << " " << r*scale << " disk" << std::endl;
			
		}
	}


	for (size_t i = 0; i < num_points; i++){
		//file<<(samples[i][0]-xmin)*scale<<" "<< (samples[i][1]-ymin)*scale<<" "<<r*scale/20.0<<" dot"<<std::endl;
#ifdef NonuniformSizingFunction
		double sizefun = SizingFunction(samples[i][0], samples[i][1], _sx, _sy, _ny, sizingFunc);
		sizefun = sqrt(sizefun);
		/*if (abs(double(_nx)*_sx - samples[i][0])<sizefun || abs(double(_ny)*_sy - samples[i][1])<sizefun ||
		samples[i][0]<sizefun || samples[i][1]<sizefun){
		continue;
		}*/
		if (samples[i][0]<XMIN || samples[i][1]<YMIN || samples[i][0]>XMAX || samples[i][1]>YMAX){
			continue;
		}
#endif

		if (use_skip_list){
			if (skip_list[i]) continue;
		}

		int length = connectivity[i][0];
		size_t j_plus = length;
		size_t ip_p = connectivity[i][length];

		for (j = 1; j <= length; j++){

			size_t ip = connectivity[i][j];
			

			if (voronoi /*&& samples[i][3]==0*/){
				if (samples[i][2] == 0){
					file << (vor[i][j * 2 - 1] - xmin)*scale << " " << (vor[i][j * 2] - ymin)*scale << " " << (vor[i][j_plus * 2 - 1] - xmin)*scale << " " << (vor[i][j_plus * 2] - ymin)*scale << " g_seg" << std::endl;
					
					j_plus = j;
				}
				else{
					file << (vor[i][j * 2 - 1] - xmin)*scale << " " << (vor[i][j * 2] - ymin)*scale << " " << (vor[i][(j + 1) * 2 - 1] - xmin)*scale << " " << (vor[i][(j + 1) * 2] - ymin)*scale << " g_seg" << std::endl;
					
				}

			}

			if (colorize){
				ip_p = FindApex(i, ip, num_points, num_points, num_points, num_points, connectivity);
				if (ip_p == num_points){
					continue;
				}
				double angles0, angles1, angles2;

				angles0 = angle_three_points_stathead(samples[ip], samples[i], samples[ip_p]);

				angles1 = angle_three_points_stathead(samples[ip_p], samples[ip], samples[i]);

				angles2 = angle_three_points_stathead(samples[i], samples[ip_p], samples[ip]);


				if (colorize && (angles0 > 90.0 + _angle_tol || angles1 > 90.0 + _angle_tol || angles2 > 90.0 + _angle_tol)){
					file << (samples[ip_p][0] - xmin)*scale << " " << (samples[ip_p][1] - ymin)*scale << " "
						<< (samples[ip][0] - xmin)*scale << " " << (samples[ip][1] - ymin)*scale << " "
						<< (samples[i][0] - xmin)*scale << " " << (samples[i][1] - ymin)*scale << " ";
					file << " red_tri" << std::endl;
				}
			}
		}

		if (voronoi /*&& samples[i][3]==0*/){
			if (samples[i][2] != 0){
				file << (vor[i][j * 2 - 1] - xmin)*scale << " " << (vor[i][j * 2] - ymin)*scale << " " << (samples[i][0] - xmin)*scale << " " << (samples[i][1] - ymin)*scale << " g_seg" << std::endl;
				file << (samples[i][0] - xmin)*scale << " " << (samples[i][1] - ymin)*scale << " " << (vor[i][(1) * 2 - 1] - xmin)*scale << " " << (vor[i][(1) * 2] - ymin)*scale << " g_seg" << std::endl;

			}
		}
	}

	if (show_connectivity){
		for (size_t i = 0; i < num_points; i++){
			if (use_skip_list){
				if (skip_list[i]){ continue; }
			}
			int length = connectivity[i][0];
			//size_t j_plus = length;
			//size_t ip_p = connectivity[i][length];

			for (int j = 1; j <= length; j++){

				size_t ip = connectivity[i][j];
				file << (samples[i][0] - xmin)*scale << " " << (samples[i][1] - ymin)*scale << " " << (samples[ip][0] - xmin)*scale << " " << (samples[ip][1] - ymin)*scale << " r_seg" << std::endl;
			}
		}
	}



	if (plot_samples_indices){
		for (size_t i = 0; i < num_points; i++){
			//if (samples[i][2] == 0) { continue; }
			if (use_skip_list){
				if (skip_list[i]) continue;
			}
			file << "/Times-Roman findfont" << std::endl;
			file << "0.08 scalefont" << std::endl;
			file << "0 0 0 setrgbcolor" << std::endl;
			file << "setfont" << std::endl;
			file << (samples[i][0] - xmin)*scale << " " << (samples[i][1] - ymin)*scale << " moveto" << std::endl;
			file << "(" << i;
			file << ") " << "show" << std::endl;
		}
	}

	if (hide_ghost_points){
		double _r = 500;
		file << (-10.0)* scale << "  " << (-10.0)*scale << "  ";
		file << (-10.0)*  scale << "  " << (ymin + 3.0*sqrt(_r)) * scale << "  ";
		file << (10.0)* scale << "  " << (ymin + 3.0*sqrt(_r))* scale << "  ";
		file << (10.0)* scale << "  " << (-10)*scale << "  ";
		file << "quad_white" << std::endl;

		file << (-10.0)* scale << "  " << (-10.0)*scale << "  ";
		file << (-10.0)*  scale << "  " << (10.0) * scale << "  ";
		file << (xmin + 3.0*sqrt(_r))* scale << "  " << (10.0)* scale << "  ";
		file << (xmin + 3.0*sqrt(_r))* scale << "  " << (-10)*scale << "  ";
		file << "quad_white" << std::endl;

		file << (-10.0)* scale << "  " << (10.0)*scale << "  ";
		file << (10.0)*  scale << "  " << (10.0) * scale << "  ";
		file << (10.0)* scale << "  " << (ymax - 3.0*sqrt(_r))* scale << "  ";
		file << (-10.0)* scale << "  " << (ymax - 3.0*sqrt(_r))*scale << "  ";
		file << "quad_white" << std::endl;

		file << (xmax - 1.0*sqrt(_r))* scale << "  " << (-10.0)*scale << "  ";
		file << (xmax - 1.0*sqrt(_r))*  scale << "  " << (10.0) * scale << "  ";
		file << (10.0)* scale << "  " << (10)* scale << "  ";
		file << (10.0)* scale << "  " << (-10)*scale << "  ";
		file << "quad_white" << std::endl;
	}
	/*for(size_t i=0;i<num_points;i++){
	if(samples[i][2]!=0){
	double min_ed_len(10E5),r_plot,dist;
	for(size_t V=1;V<=connectivity[i][0];V++){
	dist=DistPlotter(samples[i][0],samples[i][1],0,samples[connectivity[i][V]][0],samples[connectivity[i][V]][1],0);
	if(dist<min_ed_len){min_ed_len=dist;}
	}
	r_plot=sqrt(min_ed_len)/15;
	r_plot=0.001;
	if(samples[i][3]!=0){
	file<<(samples[i][0]-xmin)*scale<<" "<< (samples[i][1]-ymin)*scale<<" "<<r_plot*scale<<" purple_dot"<<std::endl;
	}else if(samples[i][2]!=0){
	file<<(samples[i][0]-xmin)*scale<<" "<< (samples[i][1]-ymin)*scale<<" "<<r_plot*scale<<" orange_dot"<<std::endl;

	}
	}
	}*/
}
void plot_grid_part(size_t n_layers, size_t ip_main, double** samples, double xo, double yo, double s, size_t* cells_i, size_t* cells_j, size_t num_cells, double r)
{
	std::fstream file("disks.ps", std::ios::app);
		

	double scale_x, scale_y, scale;

	double ss = r / sqrt(2.0);

	double xmin(samples[ip_main][0] - ss*n_layers), ymin(samples[ip_main][1] - ss*n_layers), xmax(samples[ip_main][0] + ss*n_layers), ymax(samples[ip_main][1] + ss*n_layers);

	double Lx = (xmax - xmin);
	double Ly = (ymax - ymin);

	scale_x = 6.5 / Lx;
	scale_y = 9.0 / Ly;

	if (scale_x < scale_y)
	{
		scale = scale_x;
	}
	else
	{
		scale = scale_y;
	}


	size_t ii, jj;
	for (size_t i = 0; i < num_cells; i++)
	{
		ii = cells_i[i]; jj = cells_j[i];

		file << (xo + (ii + 0) * s - xmin)* scale << "  " << (yo + (jj + 0) * s - ymin)* scale << "  ";
		file << (xo + (ii + 0) * s - xmin)* scale << "  " << (yo + (jj + 1) * s - ymin)* scale << "  ";
		file << (xo + (ii + 1) * s - xmin)* scale << "  " << (yo + (jj + 1) * s - ymin)* scale << "  ";
		file << (xo + (ii + 1) * s - xmin)* scale << "  " << (yo + (jj + 0) * s - ymin)*scale << "  ";
		file << "quad_purp" << std::endl;


	}

	file.close();
};
void plot_disk(size_t n_layers, size_t ip_main, double** samples, double xx, double yy, double r, double R)
{
	std::fstream file("disks.ps", std::ios::app);

	double scale_x, scale_y, scale;

	double ss = r / sqrt(2.0);

	double xmin(samples[ip_main][0] - ss*n_layers), ymin(samples[ip_main][1] - ss*n_layers), xmax(samples[ip_main][0] + ss*n_layers), ymax(samples[ip_main][1] + ss*n_layers);

	double Lx = (xmax - xmin);
	double Ly = (ymax - ymin);

	scale_x = 6.5 / Lx;
	scale_y = 9.0 / Ly;

	if (scale_x < scale_y)
	{
		scale = scale_x;
	}
	else
	{
		scale = scale_y;
	}

	//		
	file << (xx - xmin)*scale << " " << (yy - ymin)*scale << " " << R*scale << " blue_disk" << std::endl;
	double l = 0.15*r*cos(45.0*3.14 / 180.0);
	//
	file << (xx - l - xmin)*scale << " " << (yy - l - ymin)*scale << "	" << (xx + l - xmin)*scale << " " << (yy + l - ymin)*scale << "	r_seg" << std::endl;
	file << (xx - l - xmin)*scale << " " << (yy + l - ymin)*scale << "	" << (xx + l - xmin)*scale << " " << (yy - l - ymin)*scale << "	r_seg" << std::endl;


}


void InitiateRandNumGenerator(unsigned long x)
{

	assert(sizeof(double) >= 8);
	cc = 1.0 / 9007199254740992.0; // inverse of 2^53rd power
	int i;
	size_t qlen = indx = sizeof Q / sizeof Q[0];
	for (i = 0; i < int(qlen); i++)
		Q[i] = 0;
	int j;
	double s, t;
	unsigned long  y = 362436069;

	if (x == 0){ x = 123456789; }

	for (i = 0; i < int(qlen); i++)
	{ /* using 9th bit from Cong+Xorshift */
		s = 0.0;
		t = 1.0;
		for (j = 0; j < 52; j++)
		{
			t = 0.5 * t; /* make t=.5/2^j */
			x = 69069 * x + 123;
			y ^= (y << 13);
			y ^= (y >> 17);
			y ^= (y << 5);
			if (((x + y) >> 23) & 1)
				s = s + t; /* change bit of s, maybe */
		}	 /* end j loop */
		Q[i] = s;
	} /* end i seed loop, Now generate 10^9 RandNumGenerator's: */

}
double RandNumGenerator()
{
	double c = 0.0, zc = 0.0,	/* current CSWB and SWB `borrow` */
		zx = 5212886298506819.0 / 9007199254740992.0,	/* SWB seed1 */
		zy = 2020898595989513.0 / 9007199254740992.0;	/* SWB seed2 */

	/* Takes 14 nanosecs, Intel Q6600,2.40GHz */
	int i, j;
	double t; /* t: first temp, then next CSWB value */
	/* First get zy as next lag-2 SWB */
	t = zx - zy - zc;
	zx = zy;
	if (t < 0)
	{
		zy = t + 1.0;
		zc = cc;
	}
	else
	{
		zy = t;
		zc = 0.0;
	}

	/* Then get t as the next lag-1220 CSWB value */
	if (indx < 1220)
		t = Q[indx++];
	else
	{ /* refill Q[n] via Q[n-1220]-Q[n-1190]-c, */
		for (i = 0; i < 1220; i++)
		{
			j = (i < 30) ? i + 1190 : i - 30;
			t = Q[j] - Q[i] + c; /* Get next CSWB element */
			if (t > 0)
			{
				t = t - cc;
				c = cc;
			}
			else
			{
				t = t - cc + 1.0;
				c = 0.0;
			}
			Q[i] = t;
		}	 /* end i loop */
		indx = 1;
		t = Q[0]; /* set indx, exit 'else' with t=Q[0] */
	} /* end else segment; return t-zy mod 1 */
	return ((t < zy) ? 1.0 + (t - zy) : t - zy);
} /* end RandNumGenerator() */
double Dist(double x1, double y1, double z1, double x2, double y2, double z2)
{
	double dx, dy, dz;
	dx = x1 - x2;
	dy = y1 - y2;
	dz = z1 - z2;
	dx *= dx;
	dy *= dy;
	dz *= dz;

	return dx + dy + dz;

}
bool ConnectivityCheck(size_t node1, size_t node2)
{
	//check if node1 and node2 are connected;
	size_t V;
	for (V = 1; V <= _del[node1][0]; V++){
		if (_del[node1][V] == node2){ return true; }
	}
	return false;
}
double GetAngle360(size_t ip, size_t i1, size_t i2, size_t i3, size_t i4)
{
	double angle = 0;
	size_t ip1, ip2, d;
	ip1 = _del[ip][_del[ip][0]];
	for (d = 1; d <= _del[ip][0]; d++){
		ip2 = _del[ip][d];

		if ((ip2 == i1 || ip2 == i2 || ip2 == i3 || ip2 == i4) &&
			(ip1 == i1 || ip1 == i2 || ip1 == i3 || ip1 == i4)){
			angle += angle_three_points_stathead(_samples[ip1], _samples[ip], _samples[ip2]);
		}
		ip1 = ip2;
	}
	return angle;
}
void Partition(int& i, int& j, size_t*list)
{
	int pivot = (i + j) / 2;
	size_t tmp;
	double tmp1;
	double pivot_value = _angles[pivot];

	while (i <= j){
		while (pivot_value<_angles[j]){ j--; }
		while (pivot_value>_angles[i]){ i++; }
		if (i <= j){
			tmp = list[i + 1];
			list[i + 1] = list[j + 1];
			list[j + 1] = tmp;
			tmp1 = _angles[i];
			_angles[i] = _angles[j];
			_angles[j] = tmp1;
			j--;
			i++;
		}
	}
}
void Sort(int left, int right, size_t*list)
{
	int i, j;
	i = left;
	j = right;
	Partition(i, j, list);
	if (left< j){ Sort(left, j, list); }
	if (i<right){ Sort(i, right, list); }
}
double GetAngle(double dy, double dx)
{
	//measures angle of a vector (dx,dy) with the positive direction of x-axis
	double theta = atan2(dy, dx);
	if (theta < 0){ theta += TwoPI; }
	//	return theta*180/PI;
	return theta;
}
void SortList(double xp, double yp, size_t*list)
{
	size_t num_nodes = list[0];
	double dx, dy;
	size_t ip1;

	//double xp = _samples[ip][0];
	//double yp = _samples[ip][1];

	for (size_t i = 1; i <= num_nodes; i++)
	{
		ip1 = list[i];

		dx = _samples[ip1][0] - xp;
		dy = _samples[ip1][1] - yp;

		_angles[i - 1] = GetAngle(dy, dx);
	}
	Sort(0, num_nodes - 1, list);
}
bool IsThere(size_t entry, size_t*arr)
{
	for (size_t i = 1; i <= arr[0]; i++){
		if (arr[i] == entry){ return true; }
	}
	return false;

}
void ConnectivityBasedSortList(size_t*list, size_t**conect, bool bd)
{

	//sort the content of list around xp,yp
	//the sorting is based on the connectivity to form a connected chain
	if (bd){
		//if boundary, the first entry should be connected to one and only one on the list		
		for (size_t i = 1; i <= list[0]; i++){
			size_t ip = list[i];

			if (_samples[ip][2] == 0){ continue; }//if it's not a boundary, then go away

			size_t num_connected = 0;
			for (size_t j = 1; j <= list[0]; j++){
				if (i == j){ continue; }
				if (IsThere(list[j], conect[ip])){ num_connected++; }
			}
			if (num_connected == 1){
				//set it to the front of th list
				list[i] = list[1];
				list[1] = ip;
			}
			else if (num_connected != 2){
				std::cout << "Error (0) at ConnectivityBasedSortList(). There is something unexpected here!!!" << std::endl;
				plot_n_layers(15, ip, _samples, _del, _num_samples, 0.005, 1, 0, 1, 0, NULL, 0, 0, NULL, 0, NULL, 0, NULL, 0);
				plot_dots_n_layers_orange_list(15, ip, _samples, list, 0.005);
				//
			}
		}
	}


	for (size_t i = 1; i < list[0]; i++){
		//loop till the end -1
		int p = list[i];
		for (size_t j = i + 1; j <= list[0]; j++){
			//loop starting from next location after i till the end
			size_t q = list[j];
			if (IsThere(q, conect[p])){
				//put q infornt of p
				//move whatever infront of p to the location of q
				int tmp = list[i + 1];
				list[i + 1] = q;
				list[j] = tmp;
				break;
			}
		}
	}
}
double GetCircumcircle(double*& p1, double*& p2, double*& p3, double &xi, double &yi)
{
	//return the circumradius
	//xi,yi are the coordinates of the circumcenter

	double x1 = 0.5 * (p1[0] + p2[0]);
	double y1 = 0.5 * (p1[1] + p2[1]);

	double dx1 = p1[0] - p2[0];
	double dy1 = p1[1] - p2[1];

	double x2 = 0.5 * (p1[0] + p3[0]);
	double y2 = 0.5 * (p1[1] + p3[1]);

	double dx2 = p1[0] - p3[0];
	double dy2 = p1[1] - p3[1];

	double dx, dy;

	if (fabs(dy1) < _tol && fabs(dy2) > _tol){
		xi = x1;
		yi = (x2 - xi) * dx2 / dy2 + y2;

		dx = xi - p1[0];
		dy = yi - p1[1];

		return (dx*dx + dy*dy);
	}
	else if (fabs(dy2) < _tol && fabs(dy1)>_tol){
		xi = x2;
		yi = (x1 - xi) * dx1 / dy1 + y1;

		dx = xi - p1[0];
		dy = yi - p1[1];

		return (dx*dx + dy*dy);
	}

#ifdef debug
	else if ((fabs(dy1) < _tol  && fabs(dy2) < _tol) || (fabs(dx1) < _tol  && fabs(dx2) < _tol)){
		return(1E6);
		std::cout << " Error(0) at GetCircumcircle()... Horizontal/Vertical line .. 2D Geometry !!!!" << std::endl;
		//

	}
#endif

	double m1 = -1 * dx1 / dy1;
	double m2 = -1 * dx2 / dy2;

	xi = (y2 - y1 + m1*x1 - m2*x2) / (m1 - m2);
	yi = (x1 - xi) * dx1 / dy1 + y1;

	dx = xi - p1[0];
	dy = yi - p1[1];

	return (dx*dx + dy*dy);

}
bool EmptyCircle(double xc, double yc, double rc_2, size_t*list, size_t skip1, size_t skip2, size_t skip3, size_t skip4, size_t skip5, size_t skip6)
{
	double dist;
	for (size_t V = 1; V <= list[0]; V++){
		if (list[V] == skip1 || list[V] == skip2 || list[V] == skip3 || list[V] == skip4 || list[V] == skip5 || list[V] == skip6){ continue; }
		dist = Dist(_samples[list[V]][0], _samples[list[V]][1], 0, xc, yc, 0);
		if (dist<rc_2 - _tol*_tol){
			return false;
		}
	}
	return true;
}
size_t FindApex(size_t ip, size_t ip1, size_t skip1, size_t skip2, size_t skip3, size_t skip4)
{
	size_t i, j;

	for (i = 1; i <= _del[ip][0]; i++){
		if (_del[ip][i] == skip1 || _del[ip][i] == skip2 || _del[ip][i] == skip3 || _del[ip][i] == skip4){ continue; }
		for (j = 1; j <= _del[ip1][0]; j++){
			if (_del[ip1][j] == skip1 || _del[ip1][j] == skip2 || _del[ip][i] == skip3 || _del[ip][i] == skip4){ continue; }
			if (_del[ip][i] == _del[ip1][j]){
				return _del[ip][i];
			}
		}
	}
	return skip1;
}
bool InList(size_t ip, size_t*list){
	//check if ip is in list
	for (size_t i = 1; i <= list[0]; i++){
		if (list[i] == ip){ return true; }
	}
	return false;
}
void GetTwoBoundaryPoints(std::vector<size_t> bd, size_t ip_center, size_t&bd1, size_t&bd2)
{
	if (false){
		plot_n_layers(5, ip_center, _samples, _del, _num_samples, 0.05, 1, 0, 1, 0, NULL, 0, 0, NULL, 0, NULL, 0, NULL, 0);
	}
	size_t V, p1, p2, p1_num, p2_num;
	bool find(false), ok;
	p2 = bd[bd.size() - 1];
	for (V = 0; V<bd.size(); V++){
		p1 = bd[V];
		if (!InList(p1, _del[p2])){
			for (p1_num = 1; p1_num <= _del[ip_center][0]; p1_num++){
				if (_del[ip_center][p1_num] == p1){ break; }
			}

			for (p2_num = 1; p2_num <= _del[ip_center][0]; p2_num++){
				if (_del[ip_center][p2_num] == p2){ break; }
			}
			ok = false;

			if (p2_num == 1 || p1_num == 1 || p2_num == _del[ip_center][0] || p2_num == _del[ip_center][0]){
				if ((p2_num == 1 && p1_num == _del[ip_center][0]) ||
					(p1_num == 1 && p2_num == _del[ip_center][0]) ||
					(abs(int(p1_num) - int(p2_num)) == 1)){
					ok = true;
				}
			}
			else{
				if (abs(int(p1_num) - int(p2_num)) == 1){ ok = true; }
			}
			if (ok){
				bd1 = p1;
				bd2 = p2;
				if (!find){
					find = true;
				}
				else{
					std::cout << "Error at  GetTwoBoundaryPoints()" << std::endl;
					plot_n_layers(40, p1, _samples, _del, _num_samples, 0.005, 1, 0, 1, 0, NULL, 0, 0, NULL, 0, NULL, 0, NULL, 0);
					//
				}
			}
		}
		p2 = p1;
	}
}
bool IsCorner(size_t ip, size_t&ip1, size_t&ip2)
{
	//should send a boundary sample point
	if (_samples[ip][2] == 0){
		std::cout << "Error(0) at IsCorner(). Not an edge point." << std::endl;
		//
	}
	size_t V, d;
	std::vector<size_t> bd;

	for (V = 1; V <= _del[ip][0]; V++){
#ifdef disconnected_regions
		if (_samples[_del[ip][V]][2] != 0 && _samples[_del[ip][V]][2] == _samples[ip][2])
#else
		if (_samples[_del[ip][V]][2] != 0) 
#endif
		{
			bd.push_back(_del[ip][V]);
		}
	}
	if (bd.size()<2){
		std::cout << "Error(1) at IsCorner(). Corner discovery went wrong." << std::endl;
		plot_n_layers(40, ip, _samples, _del, _num_samples, 0.005, 1, 0, 1, 0, NULL, 0, 0, NULL, 0, NULL, 0, NULL, 0);
		plot_single_dot_n_layers(40, ip, _samples, _samples[ip][0], _samples[ip][1], 0.005);
		//
	}
	double angle1;
	for (V = 0; V<bd.size() - 1; V++){
		for (d = V + 1; d<bd.size(); d++){
			angle1 = angle_three_points_stathead(_samples[bd[V]], _samples[ip], _samples[bd[d]]);
			if (abs(angle1 - 180.0)<20){
				ip1 = bd[V];
				ip2 = bd[d];
				if (bd.size()>2){
					GetTwoBoundaryPoints(bd, ip, ip1, ip2);
				}
				return false;
			}
		}
	}
	return true;
}
bool FindExtendedSeg(size_t ip, size_t iq, size_t&ip1, size_t&ip2)
{
	//should send a boundary sample point
	if (_samples[ip][2] == 0 || _samples[iq][2] == 0){
		std::cout << "Error(0) at FindExtendedSeg(). Not an edge point." << std::endl;
		
	}
	size_t V, d;
	std::vector<size_t> bd;
	//bool find(false);
	for (V = 1; V <= _del[ip][0]; V++){
		if (_samples[_del[ip][V]][2] != 0 && _del[ip][V] != iq){
			bd.push_back(_del[ip][V]);
		}
	}
	//if(bd.size()==0){
	//	bd.push_back(iq);
	//}

	for (V = 1; V <= _del[iq][0]; V++){
		if (_samples[_del[iq][V]][2] != 0 && _del[iq][V] != ip){
			bd.push_back(_del[iq][V]);
			//find=true;
		}
	}
	//if(!find){
	//	bd.push_back(ip);
	//}

	if (bd.size()<2){
		std::cout << "Error(1) at FindExtendedSeg(). Corner discovery went wrong." << std::endl;
		
	}
	double angle1;
	for (V = 0; V<bd.size() - 1; V++){
		for (d = V + 1; d<bd.size(); d++){
			angle1 = angle_three_points_stathead(_samples[bd[V]], _samples[ip], _samples[bd[d]]);
			if (abs(angle1 - 180.0)<_smoothness_factor){
				ip1 = bd[V];
				ip2 = bd[d];
				return true;
			}
		}
	}
	return false;
}

void ReadTriangleInput(const char *node_file, const char *ele_file)
{
	//make sure that boundary markers are not supressed in the .node file 
	//numbering here always starts with zero even if the Triangle files starts with one	
	size_t tuna, V, num_ele, n1, n2, n3;
	double xmin, xmax, ymin, ymax, angle1, angle2, angle3;

	_smoothness_factor = 1.0;
	_angle_tol = 0.02;
	_smoothness_angle = 180.0 - _smoothness_factor;



	_tol = 1E-8;
	_max_active_cell = size_t(1E4);
	_min_angle_input = 360;
	_max_angle_input = 0;

	if (_first_read){
		_list = new size_t[1000];
		_circum_ex = new double*[1000];
		_circum_in = new double*[1000];
		_line_in = new double*[1000];
		_xy_org = new double*[100];
		for (V = 0; V < 100; V++){
			_xy_org[V] = new double[2];
		}
		for (V = 0; V < 1000; V++){
			_circum_ex[V] = new double[3];
			_circum_in[V] = new double[3];
			_line_in[V] = new double[8];
		}
		_active_cells_i = new size_t[_max_active_cell];
		_active_cells_j = new size_t[_max_active_cell];
		_tmp_active_cells_i = new size_t[_max_active_cell];
		_tmp_active_cells_j = new size_t[_max_active_cell];
	}


	//reading the sample coordinates	
	std::ifstream input_node, input_ele;

	input_node.open(node_file);
	input_ele.open(ele_file);

	input_node >> _num_samples >> tuna >> tuna >> tuna;

#ifdef PeriodicPoints
	if (_first_read){
		_ghost_points = new size_t*[_num_samples];
	}
	for (V = 0; V < _num_samples; V++){
		if (_first_read){
			_ghost_points[V] = new size_t[5];
		}
		_ghost_points[V][0] = 0;
	}
#endif

	size_t max_samples_size;
	max_samples_size = size_t(double(1.5)*double(_num_samples));
	if (_first_read){
		_samples = new double*[max_samples_size];
		_del = new size_t*[max_samples_size];
		_samples[_num_samples] = new double[3];
	}

	for (V = 0; V < max_samples_size; V++){
		if (_first_read){
			_samples[V] = new double[4]; //x,y,boundary_marker
		}
		_samples[V][2] = 0;
		_samples[V][3] = 0;
		if (_first_read){
			_del[V] = new size_t[_max_del];
		}
		_del[V][0] = 0;
	}

	_num_samples_input = _num_samples;
	xmin = 10E5;
	xmax = -10E5;
	ymin = 10E5;
	ymax = -10E5;


	for (V = 0; V < _num_samples; V++){
		_del[V][0] = 0;
#ifdef PeriodicPoints
		int imain;
		input_node >> tuna >> _samples[V][0] >> _samples[V][1] >> imain;
		if (_samples[V][0] > 1.0 || _samples[V][1] > 1.0 || _samples[V][0] < 0.0 || _samples[V][1] < 0.0){
			//it is a ghost
			_ghost_points[V][0] = imain;
			_ghost_points[imain][++_ghost_points[imain][0]] = V;
			_samples[V][2] = 0;
			_samples[V][3] = 0;
		}
		else{
			//it's interior and what is read is the marker
			_samples[V][2] = imain;

			if (_samples[V][0]<xmin){ xmin = _samples[V][0]; }
			if (_samples[V][0]>xmax){ xmax = _samples[V][0]; }
			if (_samples[V][1]<ymin){ ymin = _samples[V][1]; }
			if (_samples[V][1]>ymax){ ymax = _samples[V][1]; }
		}
#else				
		input_node >> tuna >> _samples[V][0] >> _samples[V][1] >> _samples[V][2];

		//if (_samples[V][2] != 0) { _samples[V][2] = 1; }
		_samples[V][3] = 0;

		if (_samples[V][0]<xmin){ xmin = _samples[V][0]; }
		if (_samples[V][0]>xmax){ xmax = _samples[V][0]; }
		if (_samples[V][1]<ymin){ ymin = _samples[V][1]; }
		if (_samples[V][1]>ymax){ ymax = _samples[V][1]; }
#endif

	}


	input_node.close();

	input_ele >> num_ele >> tuna >> tuna;



	for (V = 0; V<num_ele; V++){
		input_ele >> tuna >> n1 >> n2 >> n3;
		if (tuna != V){
			//if numbering doesn't start with zero
			n1 -= 1;
			n2 -= 1;
			n3 -= 1;
		}
#ifdef PeriodicPoints
		if (!(_samples[n1][0] > 1.0 || _samples[n1][1] > 1.0 || _samples[n1][0] < 0.0 || _samples[n1][1] < 0.0 ||
			_samples[n2][0] > 1.0 || _samples[n2][1] > 1.0 || _samples[n2][0] < 0.0 || _samples[n2][1] < 0.0 ||
			_samples[n3][0] > 1.0 || _samples[n3][1] > 1.0 || _samples[n3][0] < 0.0 || _samples[n3][1] < 0.0)){
			//get max min angle input		
			angle1 = angle_three_points_stathead(_samples[n1], _samples[n2], _samples[n3]);
			angle2 = angle_three_points_stathead(_samples[n3], _samples[n1], _samples[n2]);
			angle3 = angle_three_points_stathead(_samples[n2], _samples[n3], _samples[n1]);


			if (angle1 < _min_angle_input){ _min_angle_input = angle1; }
			if (angle2 < _min_angle_input){ _min_angle_input = angle2; }
			if (angle3 < _min_angle_input){ _min_angle_input = angle3; }

			if (angle1 > _max_angle_input){ _max_angle_input = angle1; }
			if (angle2 > _max_angle_input){ _max_angle_input = angle2; }
			if (angle3 > _max_angle_input){ _max_angle_input = angle3; }
		}

#else
		//get max min angle input		
		angle1 = angle_three_points_stathead(_samples[n1], _samples[n2], _samples[n3]);
		angle2 = angle_three_points_stathead(_samples[n3], _samples[n1], _samples[n2]);
		angle3 = angle_three_points_stathead(_samples[n2], _samples[n3], _samples[n1]);

		if (angle1<_min_angle_input){ _min_angle_input = angle1; }
		if (angle2<_min_angle_input){ _min_angle_input = angle2; }
		if (angle3<_min_angle_input){ _min_angle_input = angle3; }

		if (angle1>_max_angle_input){ _max_angle_input = angle1; }
		if (angle2>_max_angle_input){ _max_angle_input = angle2; }
		if (angle3>_max_angle_input){ _max_angle_input = angle3; }
#endif

		if (!ConnectivityCheck(n1, n2)){
			_del[n1][0]++; _del[n1][_del[n1][0]] = n2; _del[n2][0]++; _del[n2][_del[n2][0]] = n1;
		}
		if (!ConnectivityCheck(n3, n2)){
			_del[n3][0]++; _del[n3][_del[n3][0]] = n2; _del[n2][0]++; _del[n2][_del[n2][0]] = n3;
		}
		if (!ConnectivityCheck(n1, n3)){
			_del[n1][0]++; _del[n1][_del[n1][0]] = n3; _del[n3][0]++; _del[n3][_del[n3][0]] = n1;
		}
#ifdef debug
		if (_del[n1][0]>_max_del || _del[n2][0]>_max_del || _del[n3][0]>_max_del){
			std::cout << "Error(0).. increase _max_del" << std::endl;
			
		}
#endif
	}
	input_ele.close();


	//sort delaunay
	if (_first_read){
		_angles = new double[_max_del];
	}
	for (V = 0; V<_num_samples; V++){
		SortList(_samples[V][0], _samples[V][1], _del[V]);
	}

	//discover corners
	size_t ip;
	bool all_corner;

#ifdef no_boundary_sifting
	all_corner = true;
#else
	all_corner = false;
#endif // !no_boundary_sifting

	for (ip = 0; ip < _num_samples; ip++){
#ifdef PeriodicPoints
		if (_samples[ip][0] > 1.0 || _samples[ip][1] > 1.0 || _samples[ip][0] < 0.0 || _samples[ip][1] < 0.0){ continue; }
#endif
		if (_samples[ip][2] == 0){ continue; }

		if (all_corner || IsCorner(ip, n1, n2)){
			_samples[ip][3] = 1;
		}
		else{
			_samples[ip][3] = 0;
		}
	}

#ifndef PeriodicPoints
	//double scale_factor;
	////scaling inside unit square
	//for(V=0;V<_num_samples;V++){
	//	_samples[V][0]-=xmin;
	//	_samples[V][1]-=ymin;		
	//}
	//if(lx>ly){scale_factor=lx;}
	//else{scale_factor=ly;}	

	////scale_factor=1/5.8952;

	//for(V=0;V<_num_samples;V++){
	//	_samples[V][0]/=scale_factor;
	//	_samples[V][1]/=scale_factor;
	//	//_samples[V][2]/=scale_factor;
	//}
#endif

	//if (_first_read){
	//	plot_unit_box(_samples, _del, _num_samples, 0.02, 1, 0, 0, 1, 0, 0, 0, 0, NULL, 0, 0, NULL);
	//}

	_first_read = false;

}
void WriteMesh(char *node_file, char *ele_file, size_t num_samples, double min_ang, double max_ang)
{
	std::string fname_node, fname_ele;

	//directory
	fname_node = "";
	fname_ele = "";


	std::string outfilename1, outfilename2;
	std::stringstream sstm_node, sstm_ele, sstm;

	sstm_node << fname_node;
	sstm_ele << fname_ele;

	for (size_t W = 0; W<20; W++){
		if (node_file[W] != '.'){
			sstm_node << node_file[W];
			sstm << node_file[W];
			sstm_ele << node_file[W];

		}
		else{ break; }
	}

	sstm_node << "/" << sstm.str() << "_" << num_samples << "_nodes_minA_" << floor(min_ang * 100 + 0.5) / 100 << "_maxA_" << ceil(max_ang * 100 + 0.5) / 100 << ".node";
	sstm_ele << "/" << sstm.str() << "_" << num_samples << "_nodes_minA_" << floor(min_ang * 100 + 0.5) / 100 << "_maxA_" << ceil(max_ang * 100 + 0.5) / 100 << ".ele";


	std::fstream file_node(sstm_node.str().c_str(), std::ios::out);
	std::fstream file_ele(sstm_ele.str().c_str(), std::ios::out);

	file_node.precision(15);
	file_ele.precision(15);

	file_node << _num_samples << " 2 0 1 " << std::endl;
	size_t V, d, k, ip, ip1, ip2(0), num_tri(0);
	for (V = 0; V<_num_samples; V++){
		file_node << V + 1 << " " << _samples[V][0] << " " << _samples[V][1] << " " << _samples[V][2] << std::endl;
	}
	file_node.close();

	size_t** tri = new size_t*[size_t(_num_samples * 10)];
	for (V = 0; V<_num_samples * 10; V++){
		tri[V] = new size_t[3];
	}

	for (ip = 0; ip<_num_samples; ip++){
		for (V = 1; V <= _del[ip][0]; V++){
			ip1 = _del[ip][V];

			for (d = 1; d <= _del[ip1][0]; d++){
				if (_del[ip1][d] == ip){ continue; }
				for (k = 1; k <= _del[ip][0]; k++){
					if (_del[ip][k] == ip1){ continue; }
					if (_del[ip][k] == _del[ip1][d]){
						ip2 = _del[ip][k];
						if (ip<ip1 && ip<ip2 && ip1<ip2){
							tri[num_tri][0] = ip;
							tri[num_tri][1] = ip1;
							tri[num_tri][2] = ip2;
							num_tri++;
						}
					}
				}
			}
		}
	}

	file_ele << num_tri << " 3 0 " << std::endl;
	for (V = 0; V<num_tri; V++){
		file_ele << V + 1 << " " << tri[V][0] + 1 << " " << tri[V][1] + 1 << " " << tri[V][2] + 1 << std::endl;
	}
	file_ele.close();

}
bool ObtuseHead(size_t ip, size_t&ip1, size_t&ip2)
{
	size_t V, ip4, ip5;

	ip2 = _del[ip][_del[ip][0]];
	for (V = 1; V <= _del[ip][0]; V++){
		ip1 = _del[ip][V];
		ip4 = FindApex(ip1, ip2, _num_samples, _num_samples, _num_samples, _num_samples);

		if (_samples[ip1][2] != 0 && _samples[ip2][2] != 0 && !(InList(ip1, _del[ip2]))){
			ip2 = ip1;
			continue;
		}
		if (ip4 != ip){
			ip5 = FindApex(ip1, ip2, ip4, _num_samples, _num_samples, _num_samples);
			if (ip5 != ip){
				ip2 = ip1;
				continue;
			}
		}
		//ip2=FindApex(ip,ip1,_num_samples,_num_samples,_num_samples,_num_samples);
		//if(ip2==_num_samples){continue;}

		if (angle_three_points_stathead(_samples[ip1], _samples[ip], _samples[ip2]) > 90.0 + _angle_tol){
			return true;
		}
		ip2 = ip1;
	}
	return false;
}
void VerifyMeshNonObtuse(double min_ang, size_t&num_obtuse)
{
	num_obtuse = 0;
	size_t i, j, k, V, skip;
	double xc, yc, rc, angle1, angle2, angle3;
	for (i = 0; i<_num_samples; i++){
#ifdef PeriodicPoints
		if (_samples[i][0] > 1.0 || _samples[i][1] > 1.0 || _samples[i][0] < 0.0 || _samples[i][1] < 0.0){ continue; }
#endif
		for (V = 1; V <= _del[i][0]; V++){
			j = _del[i][V];
			if (V == 1){ skip = _del[i][_del[i][0]]; }
			else{ skip = _del[i][V - 1]; }
			k = FindApex(i, j, skip, _num_samples, _num_samples, _num_samples);
			if (k == skip){
				k = FindApex(i, j, _num_samples, _num_samples, _num_samples, _num_samples);
				if (k == _num_samples){ continue; }
			}
			rc = GetCircumcircle(_samples[i], _samples[j], _samples[k], xc, yc);
			if (!EmptyCircle(xc, yc, rc, _del[i], j, k, k, k, k, k) ||
				!EmptyCircle(xc, yc, rc, _del[j], i, k, k, k, k, k) ||
				!EmptyCircle(xc, yc, rc, _del[k], i, j, j, j, j, j)){
				std::cout << "Error (0) at VerifyMeshNonObtuse(). Not a delaunay triangle at (" << i << "," << j << "," << k << ")" << std::endl;
				// mark_point(_samples[i][0],_samples[i][1],_samples,NULL,_num_samples,0,0.01);
				// mark_point(_samples[j][0],_samples[j][1],_samples,NULL,_num_samples,0,0.01);
				//mark_point(_samples[k][0],_samples[k][1],_samples,NULL,_num_samples,0,0.01);
				// 				   
			}
			angle1 = angle_three_points_stathead(_samples[i], _samples[j], _samples[k]);
			angle2 = angle_three_points_stathead(_samples[k], _samples[i], _samples[j]);
			angle3 = angle_three_points_stathead(_samples[j], _samples[k], _samples[i]);
			if (angle1<min_ang - _angle_tol || angle2<min_ang - _angle_tol || angle3<min_ang - _angle_tol){
				std::cout << "Error (1) at VerifyMeshNonObtuse(). Min angle violation at (" << i << "," << j << "," << k << ")" << std::endl;
				mark_point(_samples[i][0], _samples[i][1], _samples, NULL, _num_samples, 0, 0.1);
				mark_point(_samples[j][0], _samples[j][1], _samples, NULL, _num_samples, 0, 0.1);
				mark_point(_samples[k][0], _samples[k][1], _samples, NULL, _num_samples, 0, 0.1);
				//				   
			}

			if (angle2>90 + _angle_tol){
				num_obtuse++;
				//mark_point(_samples[i][0],_samples[i][1],_samples,NULL,_num_samples,0,0.1);
				//mark_point(_samples[j][0],_samples[j][1],_samples,NULL,_num_samples,0,0.1);
				//mark_point(_samples[k][0],_samples[k][1],_samples,NULL,_num_samples,0,0.1);					 
			}
		}
	}

}

double NormalizeNumObtuseTriangles()
{
	size_t num_obtuse = 0;
	size_t total_angles = 0;
	for (size_t i = 0; i<_num_samples; i++){

#ifdef PeriodicPoints
		if (_samples[i][0] > 1.0 || _samples[i][1] > 1.0 || _samples[i][0] < 0.0 || _samples[i][1] < 0.0){ continue; }
#endif

		for (size_t V = 1; V <= _del[i][0]; V++){
			size_t j = _del[i][V];
			size_t skip;
			if (V == 1){ skip = _del[i][_del[i][0]]; }
			else{ skip = _del[i][V - 1]; }
			size_t k = FindApex(i, j, skip, _num_samples, _num_samples, _num_samples);
			if (k == skip){
				k = FindApex(i, j, _num_samples, _num_samples, _num_samples, _num_samples);
				if (k == _num_samples){ continue; }
			}
			
			double angle2 = angle_three_points_stathead(_samples[k], _samples[i], _samples[j]);		
			total_angles++;
			if (angle2>90 + _angle_tol){
				num_obtuse++;									 
			}
		}
	}
	return double(num_obtuse) / double(total_angles);
}

void AddNewNodeToDel(size_t ip1, size_t ip2, size_t ip, size_t iq)
{
	//adding ip1 to ip2 after ip and before iq
	_del[ip2][0]++;
	_del[ip2][_del[ip2][0]] = ip1;
	SortList(_samples[ip2][0], _samples[ip2][1], _del[ip2]);
	return;

	size_t V;
	bool ip_not_there(true);
	for (V = 1; V <= _del[ip2][0]; V++){
		if (_del[ip2][V] == ip){ ip_not_there = false; break; }
	}
	if (ip1 == ip2){
		std::cout << "Error(0) at AddNewNodeToDel()" << std::endl;
		
	}

	if (ip_not_there){
		_del[ip2][0]++;
		_del[ip2][_del[ip2][0]] = ip1;
		SortList(_samples[ip2][0], _samples[ip2][1], _del[ip2]);
		return;

	}


	for (V = 1; V <= _del[ip2][0]; V++){
		//sanity check if ip1 is already there
		if (_del[ip2][V] == ip1){
			return;
		}
	}

	double angle_ip2_ip, angle_ip2_ip1;
	angle_ip2_ip = GetAngle(_samples[ip][1] - _samples[ip2][1], _samples[ip][0] - _samples[ip2][0]);
	angle_ip2_ip1 = GetAngle(_samples[ip1][1] - _samples[ip2][1], _samples[ip1][0] - _samples[ip2][0]);
	//angle_ip2_iq=GetAngle(_samples[iq][1]-_samples[ip2][1],_samples[iq][0]-_samples[ip2][0]);
	if (angle_ip2_ip1<angle_ip2_ip){ ip = iq; }

	size_t i(++_del[ip2][0]);
	while (_del[ip2][i - 1] != ip){
		_del[ip2][i] = _del[ip2][i - 1];
		--i;
	}
	_del[ip2][i] = ip1;

	if (ip == iq && angle_ip2_ip1<angle_ip2_ip){
		SortList(_samples[ip2][0], _samples[ip2][1], _del[ip2]);
	}
}
void RemoveNodeFromList(size_t*list, size_t entry)
{
	//remove entry from list while preseving its sorted format
	size_t tmp, tmp1, d;
	bool find(false);
	for (d = 1; d <= list[0]; d++){
		if (list[d] == entry){ find = true; break; }
	}
	//entry is not in the list
	if (!find){ return; }
	d = list[0];
	tmp = list[d];
	while (tmp != entry){
		tmp1 = list[d - 1];
		list[d - 1] = tmp;
		tmp = tmp1;
		d--;
	}
	list[0]--;
}
void RemoveFromDel(size_t ip)
{
	size_t V, d;

	//remove ip from all the point it's connected to		
	for (V = 1; V <= _del[ip][0]; V++){
		RemoveNodeFromList(_del[_del[ip][V]], ip);
	}

	//send ip to the end of the _samples 
	//update the points connected to _num_samples-1 with ip	
	_num_samples--;
	if (ip != _num_samples){
		for (V = 0; V <= _del[_num_samples][0]; V++){
			_del[ip][V] = _del[_num_samples][V];
		}
		_samples[ip][0] = _samples[_num_samples][0];
		_samples[ip][1] = _samples[_num_samples][1];
		_samples[ip][2] = _samples[_num_samples][2];
		_samples[ip][3] = _samples[_num_samples][3];
		for (V = 1; V <= _del[_num_samples][0]; V++){
			for (d = 1; d <= _del[_del[_num_samples][V]][0]; d++){
				if (_del[_del[_num_samples][V]][d] == _num_samples){
					_del[_del[_num_samples][V]][d] = ip;
					break;
				}
			}
		}
	}



	for (V = 1; V <= _list[0]; V++){
		if (_list[V] == _num_samples){
			_list[V] = ip; break;
		}
	}
}
void MinAnglePreservMod(size_t n1, size_t n2, double min_ang, double*line, double x_mid, double y_mid)
{
	//first point on the line 
	//min_ang-=_angle_tol;
	double rot_angle = ((90 - min_ang)*PI) / 180.0, x1, y1, x2, y2;
	//line[0]=_samples[n1][0];
	//line[1]=_samples[n1][1];

	//line[2]=_samples[n2][0];
	//line[3]=_samples[n2][1];

	//first translate so as to n1 is the origin
	//rotate about n1(origin) by 90-min_angle
	//translate to the original set 
	x1 = ((_samples[n2][0] - _samples[n1][0])*cos(rot_angle) - (_samples[n2][1] - _samples[n1][1])*sin(rot_angle)) + _samples[n1][0];
	y1 = ((_samples[n2][1] - _samples[n1][1])*cos(rot_angle) + (_samples[n2][0] - _samples[n1][0])*sin(rot_angle)) + _samples[n1][1];

	//To be improved/revised
	//a naive way to check on the rotation direction
	//by doing it again with negative sign 
	//and check distanced min_point 

	//_samples[_num_samples][0]=x1;
	//_samples[_num_samples][1]=y1;
	//double dpdpdp=angle_three_points_stathead(_samples[n2],_samples[n1],_samples[_num_samples]);

	rot_angle *= -1.0;
	double x_sub = ((_samples[n2][0] - _samples[n1][0])*cos(rot_angle) - (_samples[n2][1] - _samples[n1][1])*sin(rot_angle)) + _samples[n1][0];
	double y_sub = ((_samples[n2][1] - _samples[n1][1])*cos(rot_angle) + (_samples[n2][0] - _samples[n1][0])*sin(rot_angle)) + _samples[n1][1];

	double len1(Dist(x_sub, y_sub, 0, x_mid, y_mid, 0)), len2(Dist(x1, y1, 0, x_mid, y_mid, 0));
	if (len1>len2){
		x1 = x_sub;
		y1 = y_sub;
		rot_angle *= -1.0;
	}
	line[0] = x1 - _samples[n1][0];
	line[1] = y1 - _samples[n1][1];
	line[2] = -1.0*(_samples[n1][0] * line[0] + _samples[n1][1] * line[1]);
	//translate so as to n2 is the origin
	//rotate about n2(origin) by 90-min_angle
	//translate to the original set 
	x2 = ((_samples[n1][0] - _samples[n2][0])*cos(rot_angle) - (_samples[n1][1] - _samples[n2][1])*sin(rot_angle)) + _samples[n2][0];
	y2 = ((_samples[n1][1] - _samples[n2][1])*cos(rot_angle) + (_samples[n1][0] - _samples[n2][0])*sin(rot_angle)) + _samples[n2][1];

	line[3] = x2 - _samples[n2][0];
	line[4] = y2 - _samples[n2][1];
	line[5] = -1.0*(_samples[n2][0] * line[3] + _samples[n2][1] * line[4]);

}
void AddToNeighboursList(size_t ip, size_t ip1, size_t ip2)
{
	bool there;
	for (size_t V = 1; V <= _del[ip][0]; V++){
		if (_del[ip][V] == ip1 || _del[ip][V] == ip2){ continue; }
		there = false;
		for (size_t d = 1; d <= _list[0]; d++){
			if (_list[d] == _del[ip][V]){ there = true; break; }
		}
		if (!there){
			_list[0]++;
			_list[_list[0]] = _del[ip][V];
		}
	}
}
void RefineGrid(double xo, double yo, double&s, size_t&num_active_cells)
{

	size_t icell, jcell, i, j, tmp_num_active_cells(0);
	double xmin, ymin, xmax, ymax, xfar, yfar, xnear, ynear, dx, dy;
	bool DoNotRefine = false;


	for (i = 0; i<num_active_cells; ++i){
		icell = _active_cells_i[i];
		jcell = _active_cells_j[i];
		xmin = (xo + _active_cells_i[i] * s);
		ymin = (yo + _active_cells_j[i] * s);
		xmax = (xmin + s);
		ymax = (ymin + s);

		//x_cent=(xmin+xmax)/2.0;
		//y_cent=(ymin+ymax)/2.0;

		DoNotRefine = false;

		// outside neighbouring circumcircles (exclusion regions) 
		for (j = 0; j<_num_circum_ex; j++){

			if (fabs(_circum_ex[j][0] - xmin)>fabs(_circum_ex[j][0] - xmax)){ xfar = xmin; }
			else{ xfar = xmax; }
			if (fabs(_circum_ex[j][1] - ymin)>fabs(_circum_ex[j][1] - ymax)){ yfar = ymin; }
			else{ yfar = ymax; }

			dx = _circum_ex[j][0] - xfar;
			dy = _circum_ex[j][1] - yfar;

			if ((dx*dx + dy*dy) < _circum_ex[j][2] * (1.0 + _tol)*(1.0 + _tol)){
				DoNotRefine = true;
				break;
			}
		}
		if (DoNotRefine){ continue; }







		// inside neighbouring circumcircles (inclusion regions) 
		for (j = 0; j<_num_circum_in; j++){

			if (fabs(_circum_in[j][0] - xmin)>fabs(_circum_in[j][0] - xmax)){ xnear = xmax; }
			else{ xnear = xmin; }
			if (fabs(_circum_in[j][1] - ymin)>fabs(_circum_in[j][1] - ymax)){ ynear = ymax; }
			else{ ynear = ymin; }

			dx = _circum_in[j][0] - xnear;
			dy = _circum_in[j][1] - ynear;

			if ((dx*dx + dy*dy)>_circum_in[j][2] * (1.0 - _tol)*(1.0 - _tol)){
				DoNotRefine = true;
				break;
			}
		}
		if (DoNotRefine){ continue; }







		//inside min angle preservation region
		for (j = 0; j<_num_line_in; j++){

			/*nx=_line_in[j][4]-_line_in[j][0];
			ny=_line_in[j][5]-_line_in[j][1];
			dd=-1.0*(_line_in[j][0]*nx + _line_in[j][1]*ny);
			if(x_cent*nx + y_cent*ny+ dd>-1.0*_tol){DoNotRefine=true;break;}

			nx=_line_in[j][6]-_line_in[j][2];
			ny=_line_in[j][7]-_line_in[j][3];
			dd=-1.0*(_line_in[j][2]*nx + _line_in[j][3]*ny);
			if(x_cent*nx + y_cent*ny+ dd>-1.0*_tol){DoNotRefine=true;break;}*/

			//if(x_cent*_line_in[j][0] + y_cent*_line_in[j][1]+ _line_in[j][2]>-1.0*_tol){DoNotRefine=true;break;}
			//if(x_cent*_line_in[j][3] + y_cent*_line_in[j][4]+ _line_in[j][5]>-1.0*_tol){DoNotRefine=true;break;}

			if (!((xmin*_line_in[j][0] + ymin*_line_in[j][1] + _line_in[j][2]<1.0*_tol && xmin*_line_in[j][3] + ymin*_line_in[j][4] + _line_in[j][5]<1.0*_tol) ||
				(xmax*_line_in[j][0] + ymin*_line_in[j][1] + _line_in[j][2]<1.0*_tol && xmax*_line_in[j][3] + ymin*_line_in[j][4] + _line_in[j][5]<1.0*_tol) ||
				(xmax*_line_in[j][0] + ymax*_line_in[j][1] + _line_in[j][2]<1.0*_tol && xmax*_line_in[j][3] + ymax*_line_in[j][4] + _line_in[j][5]<1.0*_tol) ||
				(xmin*_line_in[j][0] + ymax*_line_in[j][1] + _line_in[j][2]<1.0*_tol && xmin*_line_in[j][3] + ymax*_line_in[j][4] + _line_in[j][5]<1.0*_tol))){
				DoNotRefine = true; break;
			}

		}

		if (DoNotRefine){ continue; }


		icell *= 2;// move to next level of refinement 
		jcell *= 2;

		// 1st cell
		_tmp_active_cells_i[tmp_num_active_cells] = icell;
		_tmp_active_cells_j[tmp_num_active_cells] = jcell;
		tmp_num_active_cells++;

		// 2nd cell
		_tmp_active_cells_i[tmp_num_active_cells] = icell + 1;
		_tmp_active_cells_j[tmp_num_active_cells] = jcell;
		tmp_num_active_cells++;

		// 3rd cell
		_tmp_active_cells_i[tmp_num_active_cells] = icell;
		_tmp_active_cells_j[tmp_num_active_cells] = jcell + 1;
		tmp_num_active_cells++;

		// 4th cell
		_tmp_active_cells_i[tmp_num_active_cells] = icell + 1;
		_tmp_active_cells_j[tmp_num_active_cells] = jcell + 1;
		tmp_num_active_cells++;

#ifdef debug
		if (tmp_num_active_cells + 10>_max_active_cell){
			//too much refining 
			num_active_cells = _max_active_cell;//
			return;
			std::cout << "Error(0) at RefineGrid().. Grid too much refining!!!" << std::endl;
			
		}
#endif


	}

	//swapping pointers	
	num_active_cells = tmp_num_active_cells;
	for (i = 0; i<num_active_cells; ++i){
		_active_cells_i[i] = _tmp_active_cells_i[i];
		_active_cells_j[i] = _tmp_active_cells_j[i];
	}
	// move the spacing to the next level
	s *= 0.5;
}
bool CheckNewReplacement(double xp, double yp)
{

	double dx, dy, dist;
	size_t i;
	// outside neighbouring circumcircles (exclusion regions) 
	for (i = 0; i<_num_circum_ex; i++){
		dx = xp - _circum_ex[i][0];
		dy = yp - _circum_ex[i][1];
		dist = dx*dx + dy*dy;
		if (dist<_circum_ex[i][2] * (1.0 + _tol)*(1.0 + _tol)){
			return false;
		}
	}

	// inside neighbouring circumcircles (inclusion regions)	
	for (i = 0; i<_num_circum_in; i++){
		dx = xp - _circum_in[i][0];
		dy = yp - _circum_in[i][1];
		dist = dx*dx + dy*dy;
		if (dist>_circum_in[i][2] * (1.0 - _tol)*(1.0 - _tol)){
			return false;
		}
	}

	//inside min angle preservation region	
	for (i = 0; i<_num_line_in; i++){

		/*nx=_line_in[i][4]-_line_in[i][0];
		ny=_line_in[i][5]-_line_in[i][1];
		dd=-1.0*(_line_in[i][0]*nx + _line_in[i][1]*ny);
		if(xp*nx + yp*ny+ dd>-1.0*_tol){return false;}

		nx=_line_in[i][6]-_line_in[i][2];
		ny=_line_in[i][7]-_line_in[i][3];
		dd=-1.0*(_line_in[i][2]*nx + _line_in[i][3]*ny);
		if(xp*nx + yp*ny+ dd>-1.0*_tol){return false;}*/

		if (xp*_line_in[i][0] + yp*_line_in[i][1] + _line_in[i][2]>-1.0*_tol){
			return false;
		}
		if (xp*_line_in[i][3] + yp*_line_in[i][4] + _line_in[i][5]>-1.0*_tol){
			return false;
		}
	}

	return true;
}
void NoThinAngleApex(size_t n1, size_t n2, double min_angle, double x_gap, double y_gap)
{
	double x_mid((_samples[n1][0] + _samples[n2][0]) / 2.0), y_mid((_samples[n1][1] + _samples[n2][1]) / 2.0), nx(_samples[n1][1] - _samples[n2][1]), ny(-_samples[n1][0] + _samples[n2][0]);
	double len = sqrt(Dist(_samples[n1][0], _samples[n1][1], 0, _samples[n2][0], _samples[n2][1], 0));
	double prep = (len / 2.0) / (tan(0.5*min_angle*PI / 180.0));
	double x1, y1, x2, y2, factor, m;
	if (nx == 0){ m = 1E6; }
	else{ m = ny / nx; }
	factor = sqrt((prep*prep) / (1.0 + m*m));
	x1 = x_mid + factor;
	y1 = m*(x1 - x_mid) + y_mid;

	x2 = x_mid - factor;
	y2 = m*(x2 - x_mid) + y_mid;

	//i get the direction as the point the closest to ip
	//to be improved later
	if (Dist(x1, y1, 0, x_gap, y_gap, 0)<Dist(x2, y2, 0, x_gap, y_gap, 0)){
		_samples[_num_samples][0] = x1;
		_samples[_num_samples][1] = y1;
	}
	else{
		_samples[_num_samples][0] = x2;
		_samples[_num_samples][1] = y2;
	}

	_circum_in[_num_circum_in][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[_num_samples], _circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1]);
	_num_circum_in++;




}
void BoundaryEdgesCircles(size_t n1, size_t n2, size_t n1bd, size_t n2bd, double x_center, double y_center)
{
	//two circles that represents the inclusion region in which we can (re)sample and edge point and keep 
	//smooth representation of the edge 

	double x_mid((_samples[n1][0] + _samples[n2][0]) / 2.0), y_mid((_samples[n1][1] + _samples[n2][1]) / 2.0), nx(_samples[n1][1] - _samples[n2][1]), ny(-_samples[n1][0] + _samples[n2][0]);
	double len = sqrt(Dist(_samples[n1][0], _samples[n1][1], 0, _samples[n2][0], _samples[n2][1], 0));
	double prep = (len / 2.0) / (tan(0.5*_smoothness_angle*PI / 180.0));
	double x1, y1, x2, y2, factor, m;
	if (nx == 0){ m = 1E6; }
	else{ m = ny / nx; }
	factor = sqrt((prep*prep) / (1.0 + m*m));
	x1 = x_mid + factor;
	y1 = m*(x1 - x_mid) + y_mid;
	x2 = x_mid - factor;
	y2 = m*(x2 - x_mid) + y_mid;

	_samples[_num_samples][0] = x1;
	_samples[_num_samples][1] = y1;
	_circum_in[_num_circum_in][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[_num_samples], _circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1]);
	_num_circum_in++;

	_samples[_num_samples][0] = x2;
	_samples[_num_samples][1] = y2;
	_circum_in[_num_circum_in][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[_num_samples], _circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1]);
	_num_circum_in++;

	///
	double rot_angle = ((90 - _smoothness_factor)*PI) / 180.0;
	x1 = ((x_center - _samples[n1bd][0])*cos(rot_angle) - (y_center - _samples[n1bd][1])*sin(rot_angle)) + _samples[n1bd][0];
	y1 = ((y_center - _samples[n1bd][1])*cos(rot_angle) + (x_center - _samples[n1bd][0])*sin(rot_angle)) + _samples[n1bd][1];
	_line_in[_num_line_in][0] = _samples[n1bd][0] - x1;
	_line_in[_num_line_in][1] = _samples[n1bd][1] - y1;
	_line_in[_num_line_in][2] = -1.0*(_samples[n1bd][0] * _line_in[_num_line_in][0] + _samples[n1bd][1] * _line_in[_num_line_in][1]);
	rot_angle *= -1.0;
	x2 = ((x_center - _samples[n2bd][0])*cos(rot_angle) - (y_center - _samples[n2bd][1])*sin(rot_angle)) + _samples[n2bd][0];
	y2 = ((y_center - _samples[n2bd][1])*cos(rot_angle) + (x_center - _samples[n2bd][0])*sin(rot_angle)) + _samples[n2bd][1];
	_line_in[_num_line_in][3] = _samples[n2bd][0] - x2;
	_line_in[_num_line_in][4] = _samples[n2bd][1] - y2;
	_line_in[_num_line_in][5] = -1.0*(_samples[n2bd][0] * _line_in[_num_line_in][3] + _samples[n2bd][1] * _line_in[_num_line_in][4]);
	_num_line_in++;


	x2 = ((x_center - _samples[n1bd][0])*cos(rot_angle) - (y_center - _samples[n1bd][1])*sin(rot_angle)) + _samples[n1bd][0];
	y2 = ((y_center - _samples[n1bd][1])*cos(rot_angle) + (x_center - _samples[n1bd][0])*sin(rot_angle)) + _samples[n1bd][1];
	_line_in[_num_line_in][0] = _samples[n1bd][0] - x2;
	_line_in[_num_line_in][1] = _samples[n1bd][1] - y2;
	_line_in[_num_line_in][2] = -1.0*(_samples[n1bd][0] * _line_in[_num_line_in][0] + _samples[n1bd][1] * _line_in[_num_line_in][1]);
	rot_angle *= -1.0;
	x1 = ((x_center - _samples[n2bd][0])*cos(rot_angle) - (y_center - _samples[n2bd][1])*sin(rot_angle)) + _samples[n2bd][0];
	y1 = ((y_center - _samples[n2bd][1])*cos(rot_angle) + (x_center - _samples[n2bd][0])*sin(rot_angle)) + _samples[n2bd][1];
	_line_in[_num_line_in][3] = _samples[n2bd][0] - x1;
	_line_in[_num_line_in][4] = _samples[n2bd][1] - y1;
	_line_in[_num_line_in][5] = -1.0*(_samples[n2bd][0] * _line_in[_num_line_in][3] + _samples[n2bd][1] * _line_in[_num_line_in][4]);
	_num_line_in++;






	/*}else{

	double x_mid((_samples[n1][0]+_samples[n2][0])/2.0),y_mid((_samples[n1][1]+_samples[n2][1])/2.0),nx(_samples[n1][1]-_samples[n2][1]),ny(-_samples[n1][0]+_samples[n2][0]);
	double len=sqrt(Dist(_samples[n1][0],_samples[n1][1],0,_samples[n2][0],_samples[n2][1],0));
	double prep=(len/2.0)/(tan(0.5*_smoothness_angle*M_PI/180.0));
	double x1,y1,x2,y2,factor,m;
	if(nx==0){m=1E6;}
	else{m=ny/nx;}
	factor=sqrt((prep*prep)/(1.0+m*m));
	x1=x_mid+factor;
	y1=m*(x1-x_mid)+y_mid;
	x2=x_mid-factor;
	y2=m*(x2-x_mid)+y_mid;

	_samples[_num_samples][0]=x1;
	_samples[_num_samples][1]=y1;
	_circum_in[_num_circum_in][2]=GetCircumcircle(_samples[n1],_samples[n2],_samples[_num_samples],_circum_in[_num_circum_in][0],_circum_in[_num_circum_in][1]);
	_num_circum_in++;
	_samples[_num_samples][0]=x2;
	_samples[_num_samples][1]=y2;
	_circum_in[_num_circum_in][2]=GetCircumcircle(_samples[n1],_samples[n2],_samples[_num_samples],_circum_in[_num_circum_in][0],_circum_in[_num_circum_in][1]);
	_num_circum_in++;
	}*/


}

bool NonObtuseInjection(size_t ip, size_t ip1, size_t ip2, bool draw, double min_ang)
{
	//we seek to find a sample to inject so as to destroy the obtuse angle at ip (or ip1 or ip2)
	//in doing that we remove the edge ip1,ip2, but this we create a new sample that has 4-valent 
	//which is too conservative for non-obtuse
	//thus we remove another edge
	//the candidate for that are ip-ip1, ip-ip2, ip1-ip3 or ip2-ip3
	//the candidate is discareded if the gap it forms when removed is non-convex
	//we check all the candidate till we find a successful new replacement 
	double r_plot(0.04);
	size_t layer_plot(15);

	size_t ip3, n1, V, n2, ap, k, d, bd_flag;
	size_t i;
	bool no_replace;
	std::vector <size_t> candid;
	std::vector <size_t> candid_edge;

	double angle1, angle2, x_mid, y_mid;
	ip3 = FindApex(ip1, ip2, ip, ip, ip, ip);
	if (ip3 == ip){

		//cout<<"Error (0) at NonObtuseInjection()"<<endl;
		//return false;
		//
	}

	//candidate ip-ip1
	n1 = FindApex(ip, ip1, ip2, ip2, ip2, ip2);
	if (n1 != ip2){
		angle1 = angle_three_points_stathead(_samples[n1], _samples[ip], _samples[ip2]);
		angle2 = angle_three_points_stathead(_samples[n1], _samples[ip1], _samples[ip3]);
		if (angle1<180.0 - _tol && angle2<180.0 - _tol){
			candid.push_back(n1);
			candid_edge.push_back(ip);
			candid_edge.push_back(ip1);
		}
	}
	else{
		//it's a boundary edge
		candid.push_back(ip);
		candid_edge.push_back(ip);
		candid_edge.push_back(ip1);
	}

	//candidate ip-ip2
	n1 = FindApex(ip, ip2, ip1, ip1, ip1, ip1);
	if (n1 != ip1){
		angle1 = angle_three_points_stathead(_samples[n1], _samples[ip], _samples[ip1]);
		angle2 = angle_three_points_stathead(_samples[n1], _samples[ip2], _samples[ip3]);
		if (angle1<180.0 - _tol && angle2<180.0 - _tol){
			candid.push_back(n1);
			candid_edge.push_back(ip);
			candid_edge.push_back(ip2);
		}
	}
	else{
		//it's a boundary edge
		candid.push_back(ip);
		candid_edge.push_back(ip);
		candid_edge.push_back(ip2);
	}
	if (ip != ip3){
		//candidate ip1-ip3
		n1 = FindApex(ip1, ip3, ip2, ip2, ip2, ip2);
		if (n1 != ip2){
			angle1 = angle_three_points_stathead(_samples[n1], _samples[ip1], _samples[ip]);
			angle2 = angle_three_points_stathead(_samples[n1], _samples[ip3], _samples[ip2]);
			if (angle1<180.0 - _tol && angle2<180.0 - _tol){
				candid.push_back(n1);
				candid_edge.push_back(ip1);
				candid_edge.push_back(ip3);
			}
		}
		else{
			//it's a boundary edge
			candid.push_back(ip1);
			candid_edge.push_back(ip1);
			candid_edge.push_back(ip3);
		}

		//candidate ip2-ip3
		n1 = FindApex(ip2, ip3, ip1, ip1, ip1, ip1);
		if (n1 != ip1){
			angle1 = angle_three_points_stathead(_samples[n1], _samples[ip2], _samples[ip]);
			angle2 = angle_three_points_stathead(_samples[n1], _samples[ip3], _samples[ip1]);
			if (angle1<180.0 - _tol && angle2<180.0 - _tol){
				candid.push_back(n1);
				candid_edge.push_back(ip2);
				candid_edge.push_back(ip3);
			}
		}
		else{
			//it's a boundary edge
			candid.push_back(ip2);
			candid_edge.push_back(ip2);
			candid_edge.push_back(ip3);
		}
	}
	else{
		//candidate ip1-ip2
		candid.push_back(ip1);
		candid_edge.push_back(ip1);
		candid_edge.push_back(ip2);


	}

	//form the available sampling region with each candidate
	//and try to sample



	x_mid = (_samples[ip][0] + _samples[ip1][0] + _samples[ip2][0]) / 3.0;
	y_mid = (_samples[ip][1] + _samples[ip1][1] + _samples[ip2][1]) / 3.0;
	for (V = 0; V<candid.size(); V++){

		_list[1] = ip;
		_list[2] = ip1;
		_list[3] = ip2;
		if (ip == ip3){
			if (candid[V] == ip || candid[V] == ip1 || candid[V] == ip2){
				_list[0] = 3;
			}
			else{
				_list[0] = 4;
				_list[4] = candid[V];
			}

		}
		else{
			_list[4] = ip3;
			if (candid[V] == ip || candid[V] == ip1 || candid[V] == ip2 || candid[V] == ip3){
				_list[0] = 4;
			}
			else{
				_list[0] = 5;
				_list[5] = candid[V];
			}
		}
		SortList(x_mid, y_mid, _list);
#ifdef debug
		if (draw){
			plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, 0, _circum_ex, 0, NULL, 0, NULL, 1);
			plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
		}
#endif

		//find the neighbouring circumcircles (inclusion regions)
		_num_circum_in = 0;
		bd_flag = 0;
		if (_samples[ip1][2] != 0 && _samples[ip2][2] != 0 && _samples[ip1][2] == _samples[ip2][2]){
			//ip1-ip2 will be destroyed if new replacment is found
			//restricting the replacement to lie within the region that
			//keep the smooth representation of a boundary edge
			BoundaryEdgesCircles(ip1, ip2, ip1, ip2, ((_samples[ip][0] + _samples[ip1][0]) / 2.0), ((_samples[ip][1] + _samples[ip1][1]) / 2.0));
			bd_flag = size_t(_samples[ip1][2]);
		}
		if (_samples[candid_edge[2 * V]][2] != 0 && _samples[candid_edge[2 * V + 1]][2] != 0 && _list[0] == 5){
			//same thing goes for the candidate edge 
			BoundaryEdgesCircles(candid_edge[2 * V], candid_edge[2 * V + 1], candid_edge[2 * V], candid_edge[2 * V + 1], ((_samples[candid_edge[2 * V]][0] + _samples[candid_edge[2 * V + 1]][0]) / 2.0), ((_samples[candid_edge[2 * V]][1] + _samples[candid_edge[2 * V + 1]][1]) / 2.0));
			bd_flag = size_t(_samples[candid_edge[2 * V]][2]);
		}

		n2 = _list[_list[0]];
		for (i = 1; i <= _list[0]; i++){
			ap = _list[i];
			if (i == _list[0]){ n1 = _list[1]; }
			else{ n1 = _list[i + 1]; }

			NoThinAngleApex(ap, n2, min_ang, x_mid, y_mid);

			_circum_in[_num_circum_in][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[ap], _circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1]);
			if (_circum_in[_num_circum_in][2]>2){
				//if it's too large circumcircle, then there is no need to check on it
				//since it already contains the whole domain
				n2 = ap;
				continue;
			}

			if (!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[n1], n2, ap, _num_samples, _num_samples, _num_samples, _num_samples) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[n2], n1, ap, _num_samples, _num_samples, _num_samples, _num_samples) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[ap], n2, n1, _num_samples, _num_samples, _num_samples, _num_samples) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _list, n1, n2, ap, _num_samples, _num_samples, _num_samples)){
				n2 = ap; continue;
			}
			_num_circum_in++;
			n2 = ap;
		}

#ifdef debug
		if (draw){
			plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, _num_circum_in, _circum_in, 0, NULL, 0, NULL, 1);
			plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
		}
#endif

		no_replace = false;
		if (_num_circum_in>1){
			for (k = 0; k<_num_circum_in; k++){
				for (d = k + 1; d<_num_circum_in; d++){
					if (Dist(_circum_in[k][0], _circum_in[k][1], 0, _circum_in[d][0], _circum_in[d][1], 0)*(1.0 + _tol)*(1.0 + _tol)>_circum_in[k][2] + _circum_in[d][2] + 2.0*sqrt(_circum_in[k][2])*sqrt(_circum_in[d][2])){
						no_replace = true;
						break;
					}
				}
				if (no_replace){ break; }
			}
		}
		if (no_replace){
			continue;
		}

		//find the neighbouring circumcircles (exclusion regions)
		//and find inclusion area to preserve min angle (no thin triangle area)
		_num_circum_ex = 0;
		_num_line_in = 0;
		n2 = _list[_list[0]];
		for (k = 1; k <= _list[0]; k++){
			n1 = _list[k];

			if (((n1 == candid_edge[2 * V] && n2 == candid_edge[2 * V + 1]) || (n2 == candid_edge[2 * V] && n1 == candid_edge[2 * V + 1])) && (_list[0] == 5 || _list[0] == 3)){
				n2 = n1;
				continue;
			}
			ap = FindApex(n1, n2, ip, ip1, ip2, ip3);


			MinAnglePreservMod(n1, n2, min_ang, _line_in[_num_line_in], x_mid, y_mid);
			//plot_single_dot_n_layers(3,ip,_samples,_line_in[_num_line_in][4],_line_in[_num_line_in][5],0.1);
			_num_line_in++;

			//diametral circle 
			_circum_ex[_num_circum_ex][0] = (_samples[n1][0] + _samples[n2][0]) / 2.0;
			_circum_ex[_num_circum_ex][1] = (_samples[n1][1] + _samples[n2][1]) / 2.0;
			_circum_ex[_num_circum_ex][2] = Dist(_circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1], 0, _samples[n1][0], _samples[n1][1], 0);
			_num_circum_ex++;

			//non-obtuse region 
			_line_in[_num_line_in][0] = _samples[n1][0] - _samples[n2][0];
			_line_in[_num_line_in][1] = _samples[n1][1] - _samples[n2][1];
			_line_in[_num_line_in][2] = -1.0*(_samples[n1][0] * _line_in[_num_line_in][0] + _samples[n1][1] * _line_in[_num_line_in][1]);

			_line_in[_num_line_in][3] = _samples[n2][0] - _samples[n1][0];
			_line_in[_num_line_in][4] = _samples[n2][1] - _samples[n1][1];
			_line_in[_num_line_in][5] = -1.0*(_samples[n2][0] * _line_in[_num_line_in][3] + _samples[n2][1] * _line_in[_num_line_in][4]);
			_num_line_in++;

			//neighbour circumcircle
			if (ap == ip || ap == ip1){ n2 = n1; continue; }//in case there is no apex to report		
			_circum_ex[_num_circum_ex][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[ap], _circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1]);
			_num_circum_ex++;
			n2 = n1;
		}
#ifdef debug
		if (draw){
			plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, _num_circum_ex, _circum_ex, 0, NULL, 0, NULL, 1);
			plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
			for (k = 1; k <= _list[0]; k++){
				n1 = _list[k];
				ap = FindApex(n1, n2, ip, ip1, ip2, ip3);
				if (ap == ip || ap == ip1){ n2 = n1; continue; };//in case there is no apex to report
				plot_single_dot_n_layers(layer_plot, ip, _samples, _samples[ap][0], _samples[ap][1], r_plot);
				n2 = n1;
			}
		}
#endif

		//bulding the grid
		size_t num_active_cells = 16;

		// 4X4 initial cells
		_active_cells_i[0] = 0; _active_cells_j[0] = 0;
		_active_cells_i[1] = 0; _active_cells_j[1] = 1;
		_active_cells_i[2] = 0; _active_cells_j[2] = 2;
		_active_cells_i[3] = 0; _active_cells_j[3] = 3;

		_active_cells_i[4] = 1; _active_cells_j[4] = 0;
		_active_cells_i[5] = 1; _active_cells_j[5] = 1;
		_active_cells_i[6] = 1; _active_cells_j[6] = 2;
		_active_cells_i[7] = 1; _active_cells_j[7] = 3;

		_active_cells_i[8] = 2; _active_cells_j[8] = 0;
		_active_cells_i[9] = 2; _active_cells_j[9] = 1;
		_active_cells_i[10] = 2; _active_cells_j[10] = 2;
		_active_cells_i[11] = 2; _active_cells_j[11] = 3;

		_active_cells_i[12] = 3; _active_cells_j[12] = 0;
		_active_cells_i[13] = 3; _active_cells_j[13] = 1;
		_active_cells_i[14] = 3; _active_cells_j[14] = 2;
		_active_cells_i[15] = 3; _active_cells_j[15] = 3;

		double x_min_grid(10E4), y_min_grid(10E4), x_max_grid(-10E4), y_max_grid(-10E4), s, lx, ly;
		size_t idart, iactive, ii, jj;

		for (k = 1; k <= _list[0]; k++){
			if (_samples[_list[k]][0]<x_min_grid){ x_min_grid = _samples[_list[k]][0]; }
			if (_samples[_list[k]][1]<y_min_grid){ y_min_grid = _samples[_list[k]][1]; }

			if (_samples[_list[k]][0]>x_max_grid){ x_max_grid = _samples[_list[k]][0]; }
			if (_samples[_list[k]][1]>y_max_grid){ y_max_grid = _samples[_list[k]][1]; }
		}

		lx = x_max_grid - x_min_grid;
		ly = y_max_grid - y_min_grid;
		if (lx>ly){ s = lx / 4.0; }
		else{ s = ly / 4.0; }
		double x_new, y_new;


		while (true){

#ifdef debug
			if (draw){
				plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, _num_circum_ex, _circum_ex, 0, NULL, 0, NULL, 1);
				plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, _num_circum_in, _circum_in, 0, NULL, 0, NULL, 1);
				plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
				plot_grid_part(layer_plot, ip, _samples, x_min_grid, y_min_grid, s, _active_cells_i, _active_cells_j, num_active_cells, r_plot);
			}
#endif

			RefineGrid(x_min_grid, y_min_grid, s, num_active_cells);
			if (_max_active_cell == num_active_cells){
				break; //unable to find new replacement cuz it's probably lies on a line
			}
#ifdef debug
			if (draw){
				plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, _num_circum_ex, _circum_ex, 0, NULL, 0, NULL, 1);
				plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, _num_circum_in, _circum_in, 0, NULL, 0, NULL, 1);
				plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
				plot_grid_part(layer_plot, ip, _samples, x_min_grid, y_min_grid, s, _active_cells_i, _active_cells_j, num_active_cells, r_plot);
			}
#endif

			if (num_active_cells == 0){
				break;
			}

			for (idart = 0; idart<num_active_cells; idart++){
				iactive = size_t((num_active_cells - 1)*RandNumGenerator());
				ii = _active_cells_i[iactive];
				jj = _active_cells_j[iactive];

				x_new = x_min_grid + (double(ii) + RandNumGenerator())*(s);// getting the position of the dart 
				y_new = y_min_grid + (double(jj) + RandNumGenerator())*(s);


				if (CheckNewReplacement(x_new, y_new)){

					RemoveNodeFromList(_del[ip1], ip2);
					RemoveNodeFromList(_del[ip2], ip1);

					if (_list[0] != 4){
						RemoveNodeFromList(_del[candid_edge[2 * V]], candid_edge[2 * V + 1]);
						RemoveNodeFromList(_del[candid_edge[2 * V + 1]], candid_edge[2 * V]);
					}

					//update delaunay triangulation
					_samples[_num_samples][0] = x_new;
					_samples[_num_samples][1] = y_new;
					_samples[_num_samples][2] = bd_flag;

					for (k = 0; k <= _list[0]; k++){
						_del[_num_samples][k] = _list[k];
					}

					//ap=_list[_list[0]];
					for (k = 1; k <= _list[0]; k++){
						//if(k==_list[0]){n1=_list[1];}
						//else{n1=_list[k+1];}
						AddNewNodeToDel(_num_samples, _list[k], ap, n1);
						//ap=_list[k];
					}
					_num_samples++;
#ifdef debug
					if (draw){
						plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, 0, _circum_in, 0, NULL, 0, NULL, 1);
						plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
						//plot_unit_box(_samples,_del,_num_samples,0.1,0,0,0,1,0,0,0,0,NULL,0,0,NULL);
						//mark_point(_samples[ip][0],_samples[ip][1],_samples,NULL,_num_samples,0,0.03);
					}
#endif

					return true;

				}

			}

		}
	}
	return false;
}
bool NeighbourSamplesAttractor(size_t ip, size_t ip1, bool draw, double min_ang)
{
	//we seek to move the neighbour samples closer to the gap created
	//by removing ip and ip1 with an eye on the angle bounds of the 
	//triangle surround the gap
	size_t i, V, d, n2, n1, ap, i_rep, i_af, i_bf, num_new_samples;
	double x_gap(0), y_gap(0), x_best, y_best, dist, dist_best, angle;
	bool adjusted(false), find, no_replace;

	n2 = _list[_list[0]];
	for (V = 1; V <= _list[0]; V++){
		x_gap += _samples[_list[V]][0];
		y_gap += _samples[_list[V]][1];
	}
	x_gap /= _list[0];
	y_gap /= _list[0];

#ifdef debug
	if (draw){
		plot_n_layers(4, ip, _samples, _del, _num_samples, 0.05, 1, 0, 1, 0, NULL, 0, 0, _circum_ex, 0, NULL, 0, NULL, 1);
		plot_dots_n_layers_orange_list(4, ip, _samples, _list, 0.05);
	}
#endif

	for (i = 1; i <= _list[0]; i++){

		i_rep = _list[i];
		if (_samples[i_rep][2] != 0){ continue; }

		if (i == 1){ i_bf = _list[_list[0]]; }
		else{ i_bf = _list[i - 1]; }

		if (i == _list[0]){ i_af = _list[1]; }
		else{ i_af = _list[i + 1]; }

		//the inclusion region here is the circle center at (x_gap,y_gap)
		//with radius equal to the distance between (x_gap,y_gap) and i_rep
		//which ensures movment of i_rep will get it closer to the gap center
		angle = GetAngle360(i_rep, ip, ip1, i_af, i_bf);

		if (180 - angle<2.0){ continue; }


		n2 = _list[_list[0]];
		for (V = 1; V <= _list[0]; V++){
			x_gap += _samples[_list[V]][0];
			y_gap += _samples[_list[V]][1];
		}
		x_gap /= _list[0];
		y_gap /= _list[0];

		_num_circum_in = 1;

		_circum_in[0][0] = x_gap;
		_circum_in[0][1] = y_gap;
		_circum_in[0][2] = Dist(x_gap, y_gap, 0, _samples[i_rep][0], _samples[i_rep][1], 0);

		//the exclusion should maintain the delaunay edges as if we are relocating i_rep
		_num_circum_ex = 0;
		_num_line_in = 0;
		n2 = _del[i_rep][_del[i_rep][0]];

		for (V = 1; V <= _del[i_rep][0]; V++){
			n1 = _del[i_rep][V];
			if (n1 == ip || n1 == ip1 || n2 == ip || n2 == ip1){
				n2 = n1;
				continue;
			}

			ap = FindApex(n1, n2, i_rep, ip, ip1, ip1);

			NoThinAngleApex(n1, n2, min_ang, x_gap, y_gap);

			MinAnglePreservMod(n1, n2, min_ang, _line_in[_num_line_in], x_gap, y_gap);
			_num_line_in++;

			//diametral circle 
			_circum_ex[_num_circum_ex][0] = (_samples[n1][0] + _samples[n2][0]) / 2.0;
			_circum_ex[_num_circum_ex][1] = (_samples[n1][1] + _samples[n2][1]) / 2.0;
			_circum_ex[_num_circum_ex][2] = Dist(_circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1], 0, _samples[n1][0], _samples[n1][1], 0);
			_num_circum_ex++;

			//non-obtuse region 
			_line_in[_num_line_in][0] = _samples[n1][0] - _samples[n2][0];
			_line_in[_num_line_in][1] = _samples[n1][1] - _samples[n2][1];
			_line_in[_num_line_in][2] = -1.0*(_samples[n1][0] * _line_in[_num_line_in][0] + _samples[n1][1] * _line_in[_num_line_in][1]);

			_line_in[_num_line_in][3] = _samples[n2][0] - _samples[n1][0];
			_line_in[_num_line_in][4] = _samples[n2][1] - _samples[n1][1];
			_line_in[_num_line_in][5] = -1.0*(_samples[n2][0] * _line_in[_num_line_in][3] + _samples[n2][1] * _line_in[_num_line_in][4]);
			_num_line_in++;

			if (ap == i_rep){ n2 = n1; continue; }
			_circum_ex[_num_circum_ex][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[ap], _circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1]);
			_num_circum_ex++;
			n2 = n1;

		}

		n2 = _del[i_rep][_del[i_rep][0]];
		for (V = 1; V <= _del[i_rep][0]; V++){
			ap = _del[i_rep][V];
			if (V == _del[i_rep][0]){ n1 = _del[i_rep][1]; }
			else{ n1 = _del[i_rep][V + 1]; }

			if (n1 == ip || n1 == ip1 || n2 == ip || n2 == ip1 || ap == ip || ap == ip1){
				n2 = n1;
				continue;
			}

			_circum_in[_num_circum_in][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[ap], _circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1]);
			if (_circum_in[_num_circum_in][2]>2){
				//if it's too large circumcircle, then there is no need to check on it
				//since it already contains the whole domain
				n2 = ap;
				continue;
			}
			//must check if such a circle is empty
			//it could contain other existing sample point so as to it'll already never be formed
			//we check against 1-level delaunay neighbour
			if (!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[n1], n2, ap, ip, ip1, i_rep, ip1) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[n2], n1, ap, ip, ip1, i_rep, ip1) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[ap], n2, n1, ip, ip1, i_rep, ip1) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _list, n1, n2, ap, ip, ip1, i_rep)){
				n2 = ap; continue;
			}
			_num_circum_in++;
			n2 = ap;
		}

#ifdef debug
		if (draw){
			plot_n_layers(4, ip, _samples, _del, _num_samples, 0.05, 1, 0, 1, 0, NULL, 0, _num_circum_in, _circum_in, 0, NULL, 0, NULL, 1);
			plot_dots_n_layers_orange_list(4, ip, _samples, _list, 0.05);
			//mark_point(x_gap,y_gap,_samples,NULL,_num_samples,0,0.05);
			plot_disk(4, ip, _samples, x_gap, y_gap, 0.05, 0.0);
		}
#endif

		no_replace = false;
		if (_num_circum_in>1){
			for (V = 0; V<_num_circum_in; V++){
				for (d = V + 1; d<_num_circum_in; d++){
					if (Dist(_circum_in[V][0], _circum_in[V][1], 0, _circum_in[d][0], _circum_in[d][1], 0)*(1.0 + _tol)*(1.0 + _tol)>_circum_in[V][2] + _circum_in[d][2] + 2.0*sqrt(_circum_in[V][2])*sqrt(_circum_in[d][2])){
						no_replace = true;
						break;
					}
				}
				if (no_replace){ break; }
			}
		}
		if (no_replace){ continue; }


#ifdef debug
		if (draw){
			plot_n_layers(4, ip, _samples, _del, _num_samples, 0.05, 1, 0, 1, 0, NULL, 0, _num_circum_ex, _circum_ex, 0, NULL, 0, NULL, 1);
			plot_dots_n_layers_orange_list(4, ip, _samples, _list, 0.05);
			n2 = _del[i_rep][_del[i_rep][0]];
			for (V = 1; V <= _del[i_rep][0]; V++){
				n1 = _del[i_rep][V];
				if (n1 == ip || n1 == ip1 || n2 == ip || n2 == ip1){
					n2 = n1;
					continue;
				}
				ap = FindApex(n1, n2, i_rep, ip, ip1, ip1);
				if (ap == i_rep){ n2 = n1; continue; }
				plot_single_dot_n_layers(4, ip, _samples, _samples[ap][0], _samples[ap][1], 0.05);
				n2 = n1;
			}
		}
#endif

		//the exclusion region should also maintain the edges of the gap
		//the edges at are not connected to ip
		n2 = _list[_list[0]];
		for (V = 1; V <= _list[0]; V++){
			n1 = _list[V];
			if (n1 == i_rep || n2 == i_rep){
				n2 = n1;
				continue;
			}

			if (n1 == ip || n1 == ip1 || n2 == ip || n2 == ip1){
				n2 = n1;
				continue;
			}
			ap = FindApex(n1, n2, i_rep, ip, ip1, ip1);


			if (ap == i_rep){ n2 = n1; continue; }
			_circum_ex[_num_circum_ex][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[ap], _circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1]);
			_num_circum_ex++;

			n2 = n1;
		}

#ifdef debug
		if (draw){
			plot_n_layers(4, ip, _samples, _del, _num_samples, 0.05, 1, 0, 1, 0, NULL, 0, _num_circum_ex, _circum_ex, 0, NULL, 0, NULL, 1);
			plot_dots_n_layers_orange_list(4, ip, _samples, _list, 0.05);
			n2 = _del[i_rep][_del[i_rep][0]];
			for (V = 1; V <= _del[i_rep][0]; V++){
				n1 = _del[i_rep][V];
				if (n1 == ip || n1 == ip1 || n2 == ip || n2 == ip1){
					n2 = n1;
					continue;
				}
				ap = FindApex(n1, n2, i_rep, ip, ip1, ip1);
				if (ap == i_rep){ n2 = n1; continue; }
				plot_single_dot_n_layers(4, ip, _samples, _samples[ap][0], _samples[ap][1], 0.05);
				n2 = n1;
			}

			n2 = _list[_list[0]];
			for (V = 1; V <= _list[0]; V++){
				n1 = _list[V];
				if (n1 == i_rep || n2 == i_rep){
					n2 = n1;
					continue;
				}
				if (n1 == ip || n1 == ip1 || n2 == ip || n2 == ip1){
					n2 = n1;
					continue;
				}
				ap = FindApex(n1, n2, i_rep, ip, ip1, ip1);

				if (ap == i_rep){ n2 = n1; continue; }
				plot_single_dot_n_layers(4, ip, _samples, _samples[ap][0], _samples[ap][1], 0.05);
				n2 = n1;

			}
		}
#endif		

		//bulding the grid
		size_t num_active_cells = 16;

		// 4X4 initial cells
		_active_cells_i[0] = 0; _active_cells_j[0] = 0;
		_active_cells_i[1] = 0; _active_cells_j[1] = 1;
		_active_cells_i[2] = 0; _active_cells_j[2] = 2;
		_active_cells_i[3] = 0; _active_cells_j[3] = 3;

		_active_cells_i[4] = 1; _active_cells_j[4] = 0;
		_active_cells_i[5] = 1; _active_cells_j[5] = 1;
		_active_cells_i[6] = 1; _active_cells_j[6] = 2;
		_active_cells_i[7] = 1; _active_cells_j[7] = 3;

		_active_cells_i[8] = 2; _active_cells_j[8] = 0;
		_active_cells_i[9] = 2; _active_cells_j[9] = 1;
		_active_cells_i[10] = 2; _active_cells_j[10] = 2;
		_active_cells_i[11] = 2; _active_cells_j[11] = 3;

		_active_cells_i[12] = 3; _active_cells_j[12] = 0;
		_active_cells_i[13] = 3; _active_cells_j[13] = 1;
		_active_cells_i[14] = 3; _active_cells_j[14] = 2;
		_active_cells_i[15] = 3; _active_cells_j[15] = 3;

		double x_min_grid(10E4), y_min_grid(10E4), x_max_grid(-10E4), y_max_grid(-10E4), s, lx, ly;
		size_t idart, iactive, ii, jj, k;

		for (k = 1; k <= _del[i_rep][0]; k++){
			if (_samples[_del[i_rep][k]][0]<x_min_grid){ x_min_grid = _samples[_del[i_rep][k]][0]; }
			if (_samples[_del[i_rep][k]][1]<y_min_grid){ y_min_grid = _samples[_del[i_rep][k]][1]; }

			if (_samples[_del[i_rep][k]][0]>x_max_grid){ x_max_grid = _samples[_del[i_rep][k]][0]; }
			if (_samples[_del[i_rep][k]][1]>y_max_grid){ y_max_grid = _samples[_del[i_rep][k]][1]; }
		}

		for (k = 1; k <= _list[0]; k++){
			if (_samples[_list[k]][0]<x_min_grid){ x_min_grid = _samples[_list[k]][0]; }
			if (_samples[_list[k]][1]<y_min_grid){ y_min_grid = _samples[_list[k]][1]; }

			if (_samples[_list[k]][0]>x_max_grid){ x_max_grid = _samples[_list[k]][0]; }
			if (_samples[_list[k]][1]>y_max_grid){ y_max_grid = _samples[_list[k]][1]; }
		}

		lx = x_max_grid - x_min_grid;
		ly = y_max_grid - y_min_grid;
		if (lx>ly){ s = lx / 4.0; }
		else{ s = ly / 4.0; }
		double x_new, y_new;


		find = false;
		num_new_samples = 0;
		dist_best = 10E5;

		while (num_new_samples<10){

			RefineGrid(x_min_grid, y_min_grid, s, num_active_cells);
			if (_max_active_cell == num_active_cells){
				break; //unable to find new replacement cuz it's probably lies on a line
			}
#ifdef debug
			if (draw){
				plot_n_layers(4, ip, _samples, _del, _num_samples, 0.05, 1, 0, 1, 0, NULL, 0, _num_circum_ex, _circum_ex, 0, NULL, 0, NULL, 1);
				plot_n_layers(4, ip, _samples, _del, _num_samples, 0.05, 1, 0, 1, 0, NULL, 0, _num_circum_in, _circum_in, 0, NULL, 0, NULL, 1);
				plot_dots_n_layers_orange_list(4, ip, _samples, _list, 0.05);
				plot_grid_part(4, ip, _samples, x_min_grid, y_min_grid, s, _active_cells_i, _active_cells_j, num_active_cells, 0.05);
			}
#endif
			if (num_active_cells == 0){
				break;
			}

			for (idart = 0; idart<num_active_cells; idart++){
				iactive = size_t((num_active_cells - 1)*RandNumGenerator());
				ii = _active_cells_i[iactive];
				jj = _active_cells_j[iactive];

				x_new = x_min_grid + (double(ii) + RandNumGenerator())*(s);// getting the position of the dart 
				y_new = y_min_grid + (double(jj) + RandNumGenerator())*(s);



				if (CheckNewReplacement(x_new, y_new)){
					num_new_samples++;
					_samples[_num_samples][0] = x_new;
					_samples[_num_samples][1] = y_new;

					//must ensure that angle is not flat
					angle = GetAngle360(_num_samples, ip, ip1, i_bf, i_af);

					if (angle<177){
						dist = Dist(x_gap, y_gap, 0, x_new, y_new, 0);

						if (!find){
							find = true;
							dist_best = dist;
							x_best = x_new;
							y_best = y_new;
						}
						else{
							if (dist<dist_best){
								dist_best = dist;
								x_best = x_new;
								y_best = y_new;
							}
						}
						adjusted = true;
					}
				}
			}
		}

		if (find){
			_samples[i_rep][0] = x_best;
			_samples[i_rep][1] = y_best;

#ifdef debug
			if (draw){
				angle = GetAngle360(i_rep, ip, ip1, i_bf, i_af);
				plot_n_layers(4, ip, _samples, _del, _num_samples, 0.05, 1, 0, 1, 0, NULL, 0, 0, _circum_in, 0, NULL, 0, NULL, 1);
				plot_dots_n_layers_orange_list(4, ip, _samples, _list, 0.05);
			}
#endif	
		}

	}

	return adjusted;
}
bool NonObtuseSiftAttractor(size_t ip, size_t ip1, bool draw, double min_ang)
{
	if (_samples[ip][3] != 0 || _samples[ip1][3] != 0){ return false; }

	size_t V, d, n1, n2, ap, thrd, bd1, bd2, tmp;
	size_t i;
	double  x_mid((_samples[ip][0] + _samples[ip1][0]) / 2.0), y_mid((_samples[ip][1] + _samples[ip1][1]) / 2.0), angle1, angle2, angle3;
	bool no_replace, adjusted;
	_list[0] = 0;

	AddToNeighboursList(ip, ip1, ip1);
	AddToNeighboursList(ip1, ip, ip);
	SortList(x_mid, y_mid, _list);

	size_t bd_flag(0);
	if (_samples[ip][2] != 0){ bd_flag = size_t(_samples[ip][2]); }
	else if (_samples[ip1][2] != 0){ bd_flag = size_t(_samples[ip1][2]); }
	_num_line_in = 0;

	if (bd_flag != 0){
		if (_samples[ip][2] != 0 && _samples[ip1][2] != 0){
			//the starting segment here is the line segment connecting ip and ip1
			//since both are on the boundary
			if (!FindExtendedSeg(ip, ip1, bd1, bd2)){
				bd1 = ip;
				bd2 = ip1;
			}
		}
		else if (_samples[ip][2] != 0 && _samples[ip1][2] == 0){
			//the starting segment here is the bigger line segment the boundary point lies on
			if (IsCorner(ip, bd1, bd2)){
				//after sifting few point on the boundary a sample point that is not a corner 
				//becomes corner
				_samples[ip][3] = 1;
				return false;
			}
		}
		else if (_samples[ip1][2] != 0 && _samples[ip][2] == 0){
			if (IsCorner(ip1, bd1, bd2)){
				_samples[ip1][3] = 1;
				return false;
			}
		}

		/*//put bd1 and bd2 as the start and end of neighbours samples list (_list)
		while(_samples[_list[1]][2]!=bd_flag || _samples[_list[_list[0]]][2]!=bd_flag){
		n1=_list[_list[0]];
		for(i=1;i<=_list[0];i++){
		tmp=_list[i];
		_list[i]=n1;
		n1=tmp;
		}
		}*/
		//put bd1 and bd2 as the start and end of neighbours samples list (_list)

		size_t counter(0);
		while ((_samples[_list[1]][2] != bd_flag || _samples[_list[_list[0]]][2] != bd_flag) && (counter<_list[0])) {
			counter++;
			n1 = _list[_list[0]];
			for (i = 1; i <= _list[0]; i++) {
				tmp = _list[i];
				_list[i] = n1;
				n1 = tmp;
			}
		}
		if (_samples[_list[1]][2] != bd_flag || _samples[_list[_list[0]]][2] != bd_flag) {
			SortList(_samples[ip][0], _samples[ip][1], _list);
			size_t counter(0);
			while ((_samples[_list[1]][2] != bd_flag || _samples[_list[_list[0]]][2] != bd_flag) && (counter<_list[0])) {
				counter++;
				n1 = _list[_list[0]];
				for (i = 1; i <= _list[0]; i++) {
					tmp = _list[i];
					_list[i] = n1;
					n1 = tmp;
				}
			}
		}
		if (_samples[_list[1]][2] != bd_flag || _samples[_list[_list[0]]][2] != bd_flag) {
			SortList(_samples[ip1][0], _samples[ip1][1], _list);
			size_t counter(0);
			while ((_samples[_list[1]][2] != bd_flag || _samples[_list[_list[0]]][2] != bd_flag) && (counter<_list[0])) {
				counter++;
				n1 = _list[_list[0]];
				for (i = 1; i <= _list[0]; i++) {
					tmp = _list[i];
					_list[i] = n1;
					n1 = tmp;
				}
			}
		}
		if (_samples[_list[1]][2] != bd_flag || _samples[_list[_list[0]]][2] != bd_flag) {
			//giving up
			std::cout << "Error (0) at NonObtuseSiftAttractor(). I AM GIVING UP" << std::endl;
			
			return false;
		}
	}

	adjusted = false;
	for (i = 1; i <= _list[0]; i++){
		_xy_org[i][0] = _samples[_list[i]][0];
		_xy_org[i][1] = _samples[_list[i]][1];
	}
	if (NeighbourSamplesAttractor(ip, ip1, draw, min_ang)){
		adjusted = true;
	}



	_num_circum_in = 0;
	if (bd_flag != 0){
		//adding the two circles that will maintain the boundary edge 
		//from being distorded when ip moves
		BoundaryEdgesCircles(bd1, bd2, _list[1], _list[_list[0]], ((_samples[ip][0] + _samples[ip1][0]) / 2.0), ((_samples[ip][1] + _samples[ip1][1]) / 2.0));
	}
#ifdef debug
	if (draw){
		plot_n_layers(4, ip, _samples, _del, _num_samples, 0.05, 1, 0, 1, 0, NULL, 0, 0, _circum_ex, 0, NULL, 0, NULL, 1);
		plot_dots_n_layers_orange_list(4, ip, _samples, _list, 0.05);
	}
#endif

	//find the neighbouring circumcircles (inclusion regions)	
	if (bd_flag == 0){
		n2 = _list[_list[0]];
		for (i = 1; i <= _list[0]; i++){
			ap = _list[i];
			if (i == _list[0]){ n1 = _list[1]; }
			else{ n1 = _list[i + 1]; }


			//first get the circumcircle that prevent thin angle at the apex (new replacement)
			NoThinAngleApex(ap, n2, min_ang, x_mid, y_mid);

			_circum_in[_num_circum_in][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[ap], _circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1]);
			if (_circum_in[_num_circum_in][2]>2){
				//if it's too large circumcircle, then there is no need to check on it
				//since it already contains the whole domain
				n2 = ap;
				continue;
			}
			//must check if such a circle is empty
			//it could contain other existing sample point so as to it'll already never be formed
			//we check against 1-level delaunay neighbour
			if (!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[n1], n2, ap, ip, ip1, ip1, ip1) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[n2], n1, ap, ip, ip1, ip1, ip1) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[ap], n2, n1, ip, ip1, ip1, ip1) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _list, n1, n2, ap, ip, ip1, ip1)){
				n2 = ap; continue;
			}
			_num_circum_in++;
			n2 = ap;
		}
	}
	else{
		for (i = 1; i< _list[0] - 1; i++){
			n1 = _list[i];
			ap = _list[i + 1];
			n2 = _list[i + 2];

			//first get the circumcircle that prevent thin angle at the apex (new replacement)
			NoThinAngleApex(n1, ap, min_ang, _samples[ip][0], _samples[ip][1]);

			_circum_in[_num_circum_in][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[ap], _circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1]);
			if (_circum_in[_num_circum_in][2]>2){
				//if it's too large circumcircle, then there is no need to check on it
				//since it already contains the whole domain			
				continue;
			}
			//must check if such a circle is empty
			//it could contain other existing sample point so as to it'll already never be formed
			//we check against 1-level delaunay neighbour
			if (!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[n1], n2, ap, ip, ip, ip, ip) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[n2], n1, ap, ip, ip, ip, ip) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[ap], n2, n1, ip, ip, ip, ip) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _list, n1, n2, ap, ip, ip, ip)){
				continue;
			}
			_num_circum_in++;
		}
		NoThinAngleApex(ap, n2, min_ang, _samples[ip][0], _samples[ip][1]);
	}

#ifdef debug
	if (draw){
		plot_n_layers(4, ip, _samples, _del, _num_samples, 0.05, 1, 0, 1, 0, NULL, 0, _num_circum_in, _circum_in, 0, NULL, 0, NULL, 1);
		plot_dots_n_layers_orange_list(4, ip, _samples, _list, 0.05);
	}
#endif

	//if there is no area of intersection between the inclusion circle,
	//then for sure there is no replacement
	no_replace = false;
	if (_num_circum_in>1){
		for (V = 0; V<_num_circum_in; V++){
			for (d = V + 1; d<_num_circum_in; d++){
				if (Dist(_circum_in[V][0], _circum_in[V][1], 0, _circum_in[d][0], _circum_in[d][1], 0)*(1.0 + _tol)*(1.0 + _tol)>_circum_in[V][2] + _circum_in[d][2] + 2.0*sqrt(_circum_in[V][2])*sqrt(_circum_in[d][2])){
					no_replace = true;
					break;
				}
			}
			if (no_replace){ break; }
		}
	}
	if (no_replace){
		if (adjusted){
			for (i = 1; i <= _list[0]; i++){
				_samples[_list[i]][0] = _xy_org[i][0];
				_samples[_list[i]][1] = _xy_org[i][1];
			}
		}
		return false;
	}

	//find the neighbouring circumcircles (exclusion regions)
	//and find inclusion area to preserve min angle (no thin triangle area)
	_num_circum_ex = 0;



	if (bd_flag == 0){
		n2 = _list[_list[0]];
		for (V = 1; V <= _list[0]; V++){
			n1 = _list[V];
			if (_list[0] == 3){
				if (V == 3){ thrd = _list[1]; }
				else{ thrd = _list[V + 1]; }
				ap = FindApex(n1, n2, ip, ip1, thrd, ip1);
			}
			else{
				ap = FindApex(n1, n2, ip, ip1, ip1, ip1);
			}

			MinAnglePreservMod(n1, n2, min_ang, _line_in[_num_line_in], x_mid, y_mid);
			//plot_single_dot_n_layers(3,ip,_samples,_line_in[_num_line_in][4],_line_in[_num_line_in][5],0.1);

			_num_line_in++;

			//diametral circle 
			_circum_ex[_num_circum_ex][0] = (_samples[n1][0] + _samples[n2][0]) / 2.0;
			_circum_ex[_num_circum_ex][1] = (_samples[n1][1] + _samples[n2][1]) / 2.0;
			_circum_ex[_num_circum_ex][2] = Dist(_circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1], 0, _samples[n1][0], _samples[n1][1], 0);
			_num_circum_ex++;

			//non-obtuse region 
			_line_in[_num_line_in][0] = _samples[n1][0] - _samples[n2][0];
			_line_in[_num_line_in][1] = _samples[n1][1] - _samples[n2][1];
			_line_in[_num_line_in][2] = -1.0*(_samples[n1][0] * _line_in[_num_line_in][0] + _samples[n1][1] * _line_in[_num_line_in][1]);

			_line_in[_num_line_in][3] = _samples[n2][0] - _samples[n1][0];
			_line_in[_num_line_in][4] = _samples[n2][1] - _samples[n1][1];
			_line_in[_num_line_in][5] = -1.0*(_samples[n2][0] * _line_in[_num_line_in][3] + _samples[n2][1] * _line_in[_num_line_in][4]);
			_num_line_in++;

			//neighbour circumcircle
			if (ap == ip || ap == ip1){ n2 = n1; continue; }//in case there is no apex to report		
			_circum_ex[_num_circum_ex][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[ap], _circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1]);
			_num_circum_ex++;
			n2 = n1;
		}
	}
	else{
		for (V = 1; V<_list[0]; V++){
			n1 = _list[V];
			n2 = _list[V + 1];
			if (_list[0] == 3){
				if (V == 3){ thrd = _list[1]; }
				else{ thrd = _list[V + 1]; }
				ap = FindApex(n1, n2, ip, ip, thrd, ip);
			}
			else{
				ap = FindApex(n1, n2, ip, ip, ip, ip);
			}

			MinAnglePreservMod(n1, n2, min_ang, _line_in[_num_line_in], _samples[ip][0], _samples[ip][1]);
			_num_line_in++;

			//diametral circle 
			_circum_ex[_num_circum_ex][0] = (_samples[n1][0] + _samples[n2][0]) / 2.0;
			_circum_ex[_num_circum_ex][1] = (_samples[n1][1] + _samples[n2][1]) / 2.0;
			_circum_ex[_num_circum_ex][2] = Dist(_circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1], 0, _samples[n1][0], _samples[n1][1], 0);
			_num_circum_ex++;

			//non-obtuse region 
			_line_in[_num_line_in][0] = _samples[n1][0] - _samples[n2][0];
			_line_in[_num_line_in][1] = _samples[n1][1] - _samples[n2][1];
			_line_in[_num_line_in][2] = -1.0*(_samples[n1][0] * _line_in[_num_line_in][0] + _samples[n1][1] * _line_in[_num_line_in][1]);

			_line_in[_num_line_in][3] = _samples[n2][0] - _samples[n1][0];
			_line_in[_num_line_in][4] = _samples[n2][1] - _samples[n1][1];
			_line_in[_num_line_in][5] = -1.0*(_samples[n2][0] * _line_in[_num_line_in][3] + _samples[n2][1] * _line_in[_num_line_in][4]);
			_num_line_in++;

			//neighbour circumcircle
			if (ap == ip){ n2 = n1; continue; }//in case there is no apex to report		
			_circum_ex[_num_circum_ex][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[ap], _circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1]);
			_num_circum_ex++;
		}
	}

#ifdef debug
	if (draw){
		plot_n_layers(4, ip, _samples, _del, _num_samples, 0.05, 1, 0, 1, 0, NULL, 0, _num_circum_ex, _circum_ex, 0, NULL, 0, NULL, 1);
		plot_dots_n_layers_orange_list(4, ip, _samples, _list, 0.05);
		for (V = 1; V <= _list[0]; V++){
			n1 = _list[V];
			if (_list[0] == 3){
				if (V == 3){ thrd = _list[1]; }
				else{ thrd = _list[V + 1]; }
				ap = FindApex(n1, n2, ip, ip1, thrd, ip1);
			}
			else{
				ap = FindApex(n1, n2, ip, ip1, ip1, ip1);
			}
			if (ap == ip || ap == ip1){ n2 = n1; continue; };//in case there is no apex to report
			plot_single_dot_n_layers(4, ip, _samples, _samples[ap][0], _samples[ap][1], 0.05);
			n2 = n1;
		}
	}
#endif

	//bulding the grid
	size_t num_active_cells = 16;

	// 4X4 initial cells
	_active_cells_i[0] = 0; _active_cells_j[0] = 0;
	_active_cells_i[1] = 0; _active_cells_j[1] = 1;
	_active_cells_i[2] = 0; _active_cells_j[2] = 2;
	_active_cells_i[3] = 0; _active_cells_j[3] = 3;

	_active_cells_i[4] = 1; _active_cells_j[4] = 0;
	_active_cells_i[5] = 1; _active_cells_j[5] = 1;
	_active_cells_i[6] = 1; _active_cells_j[6] = 2;
	_active_cells_i[7] = 1; _active_cells_j[7] = 3;

	_active_cells_i[8] = 2; _active_cells_j[8] = 0;
	_active_cells_i[9] = 2; _active_cells_j[9] = 1;
	_active_cells_i[10] = 2; _active_cells_j[10] = 2;
	_active_cells_i[11] = 2; _active_cells_j[11] = 3;

	_active_cells_i[12] = 3; _active_cells_j[12] = 0;
	_active_cells_i[13] = 3; _active_cells_j[13] = 1;
	_active_cells_i[14] = 3; _active_cells_j[14] = 2;
	_active_cells_i[15] = 3; _active_cells_j[15] = 3;

	double x_min_grid(10E4), y_min_grid(10E4), x_max_grid(-10E4), y_max_grid(-10E4), s, lx, ly;
	size_t idart, iactive, ii, jj;

	for (V = 1; V <= _list[0]; V++){
		if (_samples[_list[V]][0]<x_min_grid){ x_min_grid = _samples[_list[V]][0]; }
		if (_samples[_list[V]][1]<y_min_grid){ y_min_grid = _samples[_list[V]][1]; }

		if (_samples[_list[V]][0]>x_max_grid){ x_max_grid = _samples[_list[V]][0]; }
		if (_samples[_list[V]][1]>y_max_grid){ y_max_grid = _samples[_list[V]][1]; }
	}

	lx = x_max_grid - x_min_grid;
	ly = y_max_grid - y_min_grid;
	if (lx>ly){ s = lx / 4.0; }
	else{ s = ly / 4.0; }

	double x_new, y_new;

	//sampling 
	bool find = false;
	double x_best, y_best, perfect_angle;
	if (bd_flag == 0){
		perfect_angle = 360.0 / double(_list[0]);
	}
	else{
		perfect_angle = 0;
		for (V = 1; V<_list[0]; V++){
			n1 = _list[V];
			n2 = _list[V + 1];
			perfect_angle += angle_three_points_stathead(_samples[n1], _samples[ip], _samples[n2]);
		}
		perfect_angle /= double(_list[0] - 1);
	}
	angle1 = 360.0;//smallest angle with the best x,y
	size_t num_new_samples(0);

	while (num_new_samples<10){

		RefineGrid(x_min_grid, y_min_grid, s, num_active_cells);
		if (_max_active_cell == num_active_cells){
			break; //unable to find new replacement cuz it's probably lies on a line
		}
#ifdef debug
		if (draw){
			plot_n_layers(4, ip, _samples, _del, _num_samples, 0.05, 1, 0, 1, 0, NULL, 0, _num_circum_ex, _circum_ex, 0, NULL, 0, NULL, 1);
			plot_n_layers(4, ip, _samples, _del, _num_samples, 0.05, 1, 0, 1, 0, NULL, 0, _num_circum_in, _circum_in, 0, NULL, 0, NULL, 1);
			plot_dots_n_layers_orange_list(4, ip, _samples, _list, 0.05);
			plot_grid_part(4, ip, _samples, x_min_grid, y_min_grid, s, _active_cells_i, _active_cells_j, num_active_cells, 0.05);
		}
#endif
		if (num_active_cells == 0){
			break;
		}

		for (idart = 0; idart<num_active_cells; idart++){
			iactive = size_t((num_active_cells - 1)*RandNumGenerator());
			ii = _active_cells_i[iactive];
			jj = _active_cells_j[iactive];

			x_new = x_min_grid + (double(ii) + RandNumGenerator())*(s);// getting the position of the dart 
			y_new = y_min_grid + (double(jj) + RandNumGenerator())*(s);

			if (CheckNewReplacement(x_new, y_new)){
				num_new_samples++;
				_samples[ip][0] = x_new;
				_samples[ip][1] = y_new;
				angle2 = 360;//smallest angle at the head with this new x,y
				if (bd_flag == 0){
					n2 = _list[_list[0]];
					for (V = 1; V <= _list[0]; V++){
						n1 = _list[V];
						angle3 = angle_three_points_stathead(_samples[n1], _samples[ip], _samples[n2]);
						if (angle3<angle2){
							angle2 = angle3;
						}
						n2 = n1;
					}
				}
				else{
					for (V = 1; V<_list[0]; V++){
						n1 = _list[V];
						n2 = _list[V + 1];
						angle3 = angle_three_points_stathead(_samples[n1], _samples[ip], _samples[n2]);
						if (angle3<angle2){
							angle2 = angle3;
						}
					}
				}

				if (!find){
					//if it's the first one
					find = true;
					x_best = x_new;
					y_best = y_new;
					angle1 = angle2;
				}
				else{
#ifdef debug
					if (perfect_angle - angle2 + _smoothness_factor<0.0 || perfect_angle - angle1 + _smoothness_factor<0){
						std::cout << "Error (0) at NonObtuseSift2()" << std::endl;
						
					}
#endif
					if (abs(perfect_angle - angle2)<abs(perfect_angle - angle1)){
						x_best = x_new;
						y_best = y_new;
						angle1 = angle2;
					}
				}
			}
		}
	}
	if (find){
		//remove ip1 from the data structure
		if (ip == _num_samples - 1){
			ip = ip1;
			ip1 = _num_samples - 1;
		}
		RemoveFromDel(ip1);

		//update delaunay triangulation
		//move ip to the new location, connect it with the neighbour nodes				
		_samples[ip][0] = x_best;
		_samples[ip][1] = y_best;
		_samples[ip][2] = bd_flag;
		for (V = 1; V <= _del[ip][0]; V++){
			RemoveNodeFromList(_del[_del[ip][V]], ip);
		}
		for (V = 0; V <= _list[0]; V++){
			_del[ip][V] = _list[V];
		}
		//ap=_list[_list[0]];
		for (V = 1; V <= _list[0]; V++){
			//if(V==_list[0]){n1=_list[1];}
			//else{n1=_list[V+1];}
			AddNewNodeToDel(ip, _list[V], ap, n1);
			//ap=_list[V];
		}

#ifdef debug
		if (draw){
			plot_n_layers(4, ip, _samples, _del, _num_samples, 0.05, 1, 0, 1, 0, NULL, 0, 0, _circum_in, 0, NULL, 0, NULL, 1);
			plot_dots_n_layers_orange_list(4, ip, _samples, _list, 0.05);
		}
#endif		
		return true;
	}
	else if (adjusted){
		//if we moved the neighbour
		//we should get them back to their original place
		for (i = 1; i <= _list[0]; i++){
			_samples[_list[i]][0] = _xy_org[i][0];
			_samples[_list[i]][1] = _xy_org[i][1];
		}
		return false;
	}

	return false;
}
bool NonObtuseSift(size_t ip, size_t ip1, bool draw, double min_ang)
{
	if (_samples[ip][3] != 0 || _samples[ip1][3] != 0){ return false; }

	size_t V, d, n1, n2, ap, thrd, bd1, bd2, tmp;
	size_t i;
	double  x_mid((_samples[ip][0] + _samples[ip1][0]) / 2.0), y_mid((_samples[ip][1] + _samples[ip1][1]) / 2.0), angle1, angle2, angle3;
	bool no_replace;
	_list[0] = 0;

	AddToNeighboursList(ip, ip1, ip1);
	AddToNeighboursList(ip1, ip, ip);
	SortList(x_mid, y_mid, _list);

	size_t bd_flag(0);
	if (_samples[ip][2] != 0){ bd_flag = size_t(_samples[ip][2]); }
	else if (_samples[ip1][2] != 0){ bd_flag = size_t(_samples[ip1][2]); }
	_num_line_in = 0;
	_num_circum_in = 0;

	if (bd_flag != 0){
		if (_samples[ip][2] != 0 && _samples[ip1][2] != 0){
			//the starting segment here is the line segment connecting ip and ip1
			//since both are on the boundary
			if (!FindExtendedSeg(ip, ip1, bd1, bd2)){
				bd1 = ip;
				bd2 = ip1;
			}
		}
		else if (_samples[ip][2] != 0 && _samples[ip1][2] == 0){
			//the starting segment here is the bigger line segment the boundary point lies on
			if (IsCorner(ip, bd1, bd2)){
				//after sifting few point on the boundary a sample point that is not a corner 
				//becomes corner
				_samples[ip][3] = 1;
				return false;
			}
		}
		else if (_samples[ip1][2] != 0 && _samples[ip][2] == 0){
			if (IsCorner(ip1, bd1, bd2)){
				_samples[ip1][3] = 1;
				return false;
			}
		}


		/*//put bd1 and bd2 as the start and end of neighbours samples list (_list)
		while(_samples[_list[1]][2]!=bd_flag || _samples[_list[_list[0]]][2]!=bd_flag){
		n1=_list[_list[0]];
		for(i=1;i<=_list[0];i++){
		tmp=_list[i];
		_list[i]=n1;
		n1=tmp;
		}
		}*/
		//put bd1 and bd2 as the start and end of neighbours samples list (_list)

		size_t counter(0);
		while ((_samples[_list[1]][2] != bd_flag || _samples[_list[_list[0]]][2] != bd_flag) && (counter<_list[0])) {
			counter++;
			n1 = _list[_list[0]];
			for (i = 1; i <= _list[0]; i++) {
				tmp = _list[i];
				_list[i] = n1;
				n1 = tmp;
			}
		}
		if (_samples[_list[1]][2] != bd_flag || _samples[_list[_list[0]]][2] != bd_flag) {
			SortList(_samples[ip][0], _samples[ip][1], _list);
			size_t counter(0);
			while ((_samples[_list[1]][2] != bd_flag || _samples[_list[_list[0]]][2] != bd_flag) && (counter<_list[0])) {
				counter++;
				n1 = _list[_list[0]];
				for (i = 1; i <= _list[0]; i++) {
					tmp = _list[i];
					_list[i] = n1;
					n1 = tmp;
				}
			}
		}
		if (_samples[_list[1]][2] != bd_flag || _samples[_list[_list[0]]][2] != bd_flag) {
			SortList(_samples[ip1][0], _samples[ip1][1], _list);
			size_t counter(0);
			while ((_samples[_list[1]][2] != bd_flag || _samples[_list[_list[0]]][2] != bd_flag) && (counter<_list[0])) {
				counter++;
				n1 = _list[_list[0]];
				for (i = 1; i <= _list[0]; i++) {
					tmp = _list[i];
					_list[i] = n1;
					n1 = tmp;
				}
			}
		}
		if (_samples[_list[1]][2] != bd_flag || _samples[_list[_list[0]]][2] != bd_flag) {
			//giving up
			std::cout << "Error (0) at NonObtuseSift2(). I AM GIVING UP" << std::endl;
			
			return false;
		}
		//adding the two circles that will maintain the boundary edge 
		//from being distorded when ip moves
		BoundaryEdgesCircles(bd1, bd2, _list[1], _list[_list[0]], ((_samples[ip][0] + _samples[ip1][0]) / 2.0), ((_samples[ip][1] + _samples[ip1][1]) / 2.0));
	}

#ifdef debug
	if (draw){
		plot_n_layers(4, ip, _samples, _del, _num_samples, 0.05, 1, 0, 1, 0, NULL, 0, 0, _circum_ex, 0, NULL, 0, NULL, 1);
		plot_dots_n_layers_orange_list(4, ip, _samples, _list, 0.05);
	}
#endif

	//find the neighbouring circumcircles (inclusion regions)	
	if (bd_flag == 0){
		n2 = _list[_list[0]];
		for (i = 1; i <= _list[0]; i++){
			ap = _list[i];
			if (i == _list[0]){ n1 = _list[1]; }
			else{ n1 = _list[i + 1]; }


			//first get the circumcircle that prevent thin angle at the apex (new replacement)
			NoThinAngleApex(ap, n2, min_ang, x_mid, y_mid);

			_circum_in[_num_circum_in][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[ap], _circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1]);
			if (_circum_in[_num_circum_in][2]>2){
				//if it's too large circumcircle, then there is no need to check on it
				//since it already contains the whole domain
				n2 = ap;
				continue;
			}
			//must check if such a circle is empty
			//it could contain other existing sample point so as to it'll already never be formed
			//we check against 1-level delaunay neighbour
			if (!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[n1], n2, ap, ip, ip1, ip1, ip1) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[n2], n1, ap, ip, ip1, ip1, ip1) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[ap], n2, n1, ip, ip1, ip1, ip1) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _list, n1, n2, ap, ip, ip1, ip1)){
				n2 = ap; continue;
			}
			_num_circum_in++;
			n2 = ap;
		}
	}
	else{
		for (i = 1; i< _list[0] - 1; i++){
			n1 = _list[i];
			ap = _list[i + 1];
			n2 = _list[i + 2];

			//first get the circumcircle that prevent thin angle at the apex (new replacement)
			NoThinAngleApex(n1, ap, min_ang, _samples[ip][0], _samples[ip][1]);

			_circum_in[_num_circum_in][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[ap], _circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1]);
			if (_circum_in[_num_circum_in][2]>2){
				//if it's too large circumcircle, then there is no need to check on it
				//since it already contains the whole domain			
				continue;
			}
			//must check if such a circle is empty
			//it could contain other existing sample point so as to it'll already never be formed
			//we check against 1-level delaunay neighbour
			if (!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[n1], n2, ap, ip, ip1, ip1, ip1) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[n2], n1, ap, ip, ip1, ip1, ip1) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[ap], n2, n1, ip, ip1, ip1, ip1) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _list, n1, n2, ap, ip, ip1, ip1)){
				continue;
			}
			_num_circum_in++;
		}
		NoThinAngleApex(ap, n2, min_ang, _samples[ip][0], _samples[ip][1]);
	}

#ifdef debug
	if (draw){
		plot_n_layers(4, ip, _samples, _del, _num_samples, 0.05, 1, 0, 1, 0, NULL, 0, _num_circum_in, _circum_in, 0, NULL, 0, NULL, 1);
		plot_dots_n_layers_orange_list(4, ip, _samples, _list, 0.05);
	}
#endif

	//if there is no area of intersection between the inclusion circle,
	//then for sure there is no replacement
	no_replace = false;
	if (_num_circum_in>1){
		for (V = 0; V<_num_circum_in; V++){
			for (d = V + 1; d<_num_circum_in; d++){
				if (Dist(_circum_in[V][0], _circum_in[V][1], 0, _circum_in[d][0], _circum_in[d][1], 0)*(1.0 + _tol)*(1.0 + _tol)>_circum_in[V][2] + _circum_in[d][2] + 2.0*sqrt(_circum_in[V][2])*sqrt(_circum_in[d][2])){
					no_replace = true;
					break;
				}
			}
			if (no_replace){ break; }
		}
	}
	if (no_replace){
		return false;
	}

	//find the neighbouring circumcircles (exclusion regions)
	//and find inclusion area to preserve min angle (no thin triangle area)
	_num_circum_ex = 0;



	if (bd_flag == 0){
		n2 = _list[_list[0]];
		for (V = 1; V <= _list[0]; V++){
			n1 = _list[V];
			if (_list[0] == 3){
				if (V == 3){ thrd = _list[1]; }
				else{ thrd = _list[V + 1]; }
				ap = FindApex(n1, n2, ip, ip1, thrd, ip1);
			}
			else{
				ap = FindApex(n1, n2, ip, ip1, ip1, ip1);
			}

			MinAnglePreservMod(n1, n2, min_ang, _line_in[_num_line_in], x_mid, y_mid);
			//plot_single_dot_n_layers(3,ip,_samples,_line_in[_num_line_in][4],_line_in[_num_line_in][5],0.1);			
			_num_line_in++;

			//diametral circle 
			_circum_ex[_num_circum_ex][0] = (_samples[n1][0] + _samples[n2][0]) / 2.0;
			_circum_ex[_num_circum_ex][1] = (_samples[n1][1] + _samples[n2][1]) / 2.0;
			_circum_ex[_num_circum_ex][2] = Dist(_circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1], 0, _samples[n1][0], _samples[n1][1], 0);
			_num_circum_ex++;

			//non-obtuse region 
			_line_in[_num_line_in][0] = _samples[n1][0] - _samples[n2][0];
			_line_in[_num_line_in][1] = _samples[n1][1] - _samples[n2][1];
			_line_in[_num_line_in][2] = -1.0*(_samples[n1][0] * _line_in[_num_line_in][0] + _samples[n1][1] * _line_in[_num_line_in][1]);

			_line_in[_num_line_in][3] = _samples[n2][0] - _samples[n1][0];
			_line_in[_num_line_in][4] = _samples[n2][1] - _samples[n1][1];
			_line_in[_num_line_in][5] = -1.0*(_samples[n2][0] * _line_in[_num_line_in][3] + _samples[n2][1] * _line_in[_num_line_in][4]);
			_num_line_in++;

			//neighbour circumcircle
			if (ap == ip || ap == ip1){ n2 = n1; continue; }//in case there is no apex to report		
			_circum_ex[_num_circum_ex][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[ap], _circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1]);
			_num_circum_ex++;
			n2 = n1;
		}
	}
	else{
		for (V = 1; V<_list[0]; V++){
			n1 = _list[V];
			n2 = _list[V + 1];
			if (_list[0] == 3){
				if (V == 3){ thrd = _list[1]; }
				else{ thrd = _list[V + 1]; }
				ap = FindApex(n1, n2, ip, ip, thrd, ip);
			}
			else{
				ap = FindApex(n1, n2, ip, ip, ip, ip);
			}

			MinAnglePreservMod(n1, n2, min_ang, _line_in[_num_line_in], _samples[ip][0], _samples[ip][1]);
			_num_line_in++;

			//diametral circle 
			_circum_ex[_num_circum_ex][0] = (_samples[n1][0] + _samples[n2][0]) / 2.0;
			_circum_ex[_num_circum_ex][1] = (_samples[n1][1] + _samples[n2][1]) / 2.0;
			_circum_ex[_num_circum_ex][2] = Dist(_circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1], 0, _samples[n1][0], _samples[n1][1], 0);
			_num_circum_ex++;

			//non-obtuse region 
			_line_in[_num_line_in][0] = _samples[n1][0] - _samples[n2][0];
			_line_in[_num_line_in][1] = _samples[n1][1] - _samples[n2][1];
			_line_in[_num_line_in][2] = -1.0*(_samples[n1][0] * _line_in[_num_line_in][0] + _samples[n1][1] * _line_in[_num_line_in][1]);

			_line_in[_num_line_in][3] = _samples[n2][0] - _samples[n1][0];
			_line_in[_num_line_in][4] = _samples[n2][1] - _samples[n1][1];
			_line_in[_num_line_in][5] = -1.0*(_samples[n2][0] * _line_in[_num_line_in][3] + _samples[n2][1] * _line_in[_num_line_in][4]);
			_num_line_in++;

			//neighbour circumcircle
			if (ap == ip){ n2 = n1; continue; }//in case there is no apex to report		
			_circum_ex[_num_circum_ex][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[ap], _circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1]);
			_num_circum_ex++;
		}
	}

#ifdef debug
	if (draw){
		plot_n_layers(4, ip, _samples, _del, _num_samples, 0.05, 1, 0, 1, 0, NULL, 0, _num_circum_ex, _circum_ex, 0, NULL, 0, NULL, 1);
		plot_dots_n_layers_orange_list(4, ip, _samples, _list, 0.05);
		for (V = 1; V <= _list[0]; V++){
			n1 = _list[V];
			if (_list[0] == 3){
				if (V == 3){ thrd = _list[1]; }
				else{ thrd = _list[V + 1]; }
				ap = FindApex(n1, n2, ip, ip1, thrd, ip1);
			}
			else{
				ap = FindApex(n1, n2, ip, ip1, ip1, ip1);
			}
			if (ap == ip || ap == ip1){ n2 = n1; continue; };//in case there is no apex to report
			plot_single_dot_n_layers(4, ip, _samples, _samples[ap][0], _samples[ap][1], 0.05);
			n2 = n1;
		}
	}
#endif

	//bulding the grid
	size_t num_active_cells = 16;

	// 4X4 initial cells
	_active_cells_i[0] = 0; _active_cells_j[0] = 0;
	_active_cells_i[1] = 0; _active_cells_j[1] = 1;
	_active_cells_i[2] = 0; _active_cells_j[2] = 2;
	_active_cells_i[3] = 0; _active_cells_j[3] = 3;

	_active_cells_i[4] = 1; _active_cells_j[4] = 0;
	_active_cells_i[5] = 1; _active_cells_j[5] = 1;
	_active_cells_i[6] = 1; _active_cells_j[6] = 2;
	_active_cells_i[7] = 1; _active_cells_j[7] = 3;

	_active_cells_i[8] = 2; _active_cells_j[8] = 0;
	_active_cells_i[9] = 2; _active_cells_j[9] = 1;
	_active_cells_i[10] = 2; _active_cells_j[10] = 2;
	_active_cells_i[11] = 2; _active_cells_j[11] = 3;

	_active_cells_i[12] = 3; _active_cells_j[12] = 0;
	_active_cells_i[13] = 3; _active_cells_j[13] = 1;
	_active_cells_i[14] = 3; _active_cells_j[14] = 2;
	_active_cells_i[15] = 3; _active_cells_j[15] = 3;

	double x_min_grid(10E4), y_min_grid(10E4), x_max_grid(-10E4), y_max_grid(-10E4), s, lx, ly;
	size_t idart, iactive, ii, jj;

	for (V = 1; V <= _list[0]; V++){
		if (_samples[_list[V]][0]<x_min_grid){ x_min_grid = _samples[_list[V]][0]; }
		if (_samples[_list[V]][1]<y_min_grid){ y_min_grid = _samples[_list[V]][1]; }

		if (_samples[_list[V]][0]>x_max_grid){ x_max_grid = _samples[_list[V]][0]; }
		if (_samples[_list[V]][1]>y_max_grid){ y_max_grid = _samples[_list[V]][1]; }
	}

	lx = x_max_grid - x_min_grid;
	ly = y_max_grid - y_min_grid;
	if (lx>ly){ s = lx / 4.0; }
	else{ s = ly / 4.0; }

	double x_new, y_new;

	//sampling 
	bool find = false;
	double x_best, y_best, perfect_angle;
	if (bd_flag == 0){
		perfect_angle = 360.0 / double(_list[0]);
	}
	else{
		perfect_angle = 0;
		for (V = 1; V<_list[0]; V++){
			n1 = _list[V];
			n2 = _list[V + 1];
			perfect_angle += angle_three_points_stathead(_samples[n1], _samples[ip], _samples[n2]);
		}
		perfect_angle /= double(_list[0] - 1);
	}
	angle1 = 360.0;//smallest angle with the best x,y
	size_t num_new_samples(0);

	while (num_new_samples<10){

		RefineGrid(x_min_grid, y_min_grid, s, num_active_cells);
		if (_max_active_cell == num_active_cells){
			return false; //unable to find new replacement cuz it's probably lies on a line
		}
#ifdef debug
		if (draw){
			plot_n_layers(4, ip, _samples, _del, _num_samples, 0.05, 1, 0, 1, 0, NULL, 0, _num_circum_ex, _circum_ex, 0, NULL, 0, NULL, 1);
			plot_n_layers(4, ip, _samples, _del, _num_samples, 0.05, 1, 0, 1, 0, NULL, 0, _num_circum_in, _circum_in, 0, NULL, 0, NULL, 1);
			plot_dots_n_layers_orange_list(4, ip, _samples, _list, 0.05);
			plot_grid_part(4, ip, _samples, x_min_grid, y_min_grid, s, _active_cells_i, _active_cells_j, num_active_cells, 0.05);
		}
#endif
		if (num_active_cells == 0){
			return false;
		}

		for (idart = 0; idart<10 * num_active_cells; idart++){
			iactive = size_t((num_active_cells - 1)*RandNumGenerator());
			ii = _active_cells_i[iactive];
			jj = _active_cells_j[iactive];

			x_new = x_min_grid + (double(ii) + RandNumGenerator())*(s);// getting the position of the dart 
			y_new = y_min_grid + (double(jj) + RandNumGenerator())*(s);

			//plot_single_dot_n_layers(4, ip, _samples, x_new, y_new, 0.05);

			if (CheckNewReplacement(x_new, y_new)){
				num_new_samples++;
				_samples[ip][0] = x_new;
				_samples[ip][1] = y_new;
				angle2 = 360;//smallest angle at the head with this new x,y
				if (bd_flag == 0){
					n2 = _list[_list[0]];
					for (V = 1; V <= _list[0]; V++){
						n1 = _list[V];
						angle3 = angle_three_points_stathead(_samples[n1], _samples[ip], _samples[n2]);
						if (angle3<angle2){
							angle2 = angle3;
						}
						n2 = n1;
					}
				}
				else{
					for (V = 1; V<_list[0]; V++){
						n1 = _list[V];
						n2 = _list[V + 1];
						angle3 = angle_three_points_stathead(_samples[n1], _samples[ip], _samples[n2]);
						if (angle3<angle2){
							angle2 = angle3;
						}
					}
				}

				if (!find){
					//if it's the first one
					find = true;
					x_best = x_new;
					y_best = y_new;
					angle1 = angle2;
				}
				else{
#ifdef debug

					if (perfect_angle - angle2 + _smoothness_factor<0.0 || perfect_angle - angle1 + _smoothness_factor<0){
						std::cout << "Error (0) at NonObtuseSift2()" << std::endl;
						
					}
#endif
					if (abs(perfect_angle - angle2)<abs(perfect_angle - angle1)){
						x_best = x_new;
						y_best = y_new;
						angle1 = angle2;
					}
				}
			}
		}
	}
	if (find){
		//remove ip1 from the data structure
		if (ip == _num_samples - 1){
			ip = ip1;
			ip1 = _num_samples - 1;
		}
		RemoveFromDel(ip1);

		//update delaunay triangulation
		//move ip to the new location, connect it with the neighbour nodes				
		_samples[ip][0] = x_best;
		_samples[ip][1] = y_best;
		_samples[ip][2] = bd_flag;
		for (V = 1; V <= _del[ip][0]; V++){
			RemoveNodeFromList(_del[_del[ip][V]], ip);
		}
		for (V = 0; V <= _list[0]; V++){
			_del[ip][V] = _list[V];
		}
		//ap=_list[_list[0]];
		for (V = 1; V <= _list[0]; V++){
			//if(V==_list[0]){n1=_list[1];}
			//else{n1=_list[V+1];}
			AddNewNodeToDel(ip, _list[V], ap, n1);
			//ap=_list[V];
		}

#ifdef debug
		if (draw){
			plot_n_layers(4, ip, _samples, _del, _num_samples, 0.05, 1, 0, 1, 0, NULL, 0, 0, _circum_in, 0, NULL, 0, NULL, 1);
			plot_dots_n_layers_orange_list(4, ip, _samples, _list, 0.05);
		}
#endif		
		return true;
	}

	return false;
}
bool NeighbourSamplesRepeller(size_t ip, size_t ip1, size_t ip2, size_t ip3, bool draw, double min_ang, size_t layer_plot, double r_plot)
{
	//spicify min edge lenght in the gap
	size_t i, V, n2, n1, i_rep, i_bf, i_af, ap, num_new_samples;
	double x_gap(0), y_gap(0), min_local_angle, x_best, y_best, angle_best, angle;
	bool adjusted(false), find;

	n2 = _list[_list[0]];
	for (V = 1; V <= _list[0]; V++){
		x_gap += _samples[_list[V]][0];
		y_gap += _samples[_list[V]][1];
	}
	x_gap /= _list[0];
	y_gap /= _list[0];


#ifdef debug
	if (draw){
		plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, 0, _circum_ex, 0, NULL, 0, NULL, 0);
		plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
	}
#endif
	size_t i_intr;
	for (i = 1; i <= _list[0]; i++){

		i_rep = _list[i];
		if (_samples[i_rep][2] != 0){ continue; }

		//if(i_rep==min_edg1 || i_rep==min_edg2){continue;}


		if (i == 1){ i_bf = _list[_list[0]]; }
		else{ i_bf = _list[i - 1]; }

		if (i == _list[0]){ i_af = _list[1]; }
		else{ i_af = _list[i + 1]; }

		//first exclusion region is a circle inside which the angle will increase 
		if (i_rep == ip || i_rep == ip1 || i_rep == ip2 || i_rep == ip3){
			if (i_rep == ip){ i_intr = ip1; }
			else if (i_rep == ip1){ i_intr = ip; }
			else if (i_rep == ip2){ i_intr = ip3; }
			else if (i_rep == ip3){ i_intr = ip2; }

			min_local_angle = angle_three_points_stathead(_samples[i_bf], _samples[i_rep], _samples[i_intr]);
			min_local_angle += angle_three_points_stathead(_samples[i_intr], _samples[i_rep], _samples[i_af]);

		}
		else{
			min_local_angle = angle_three_points_stathead(_samples[i_bf], _samples[i_rep], _samples[i_af]);
		}


		//we only seek to decrease the angle of corner with large angles 
		//because decreasing the small angles, may increase the void edges length
		//which decreases the area available for sampling
		//so if it's alreay small, we skip it

		if ((min_local_angle<90 && _list[0] == 4) || (min_local_angle<110 && _list[0] == 5)){ continue; }

		//to maure it moves away from the void 	

		_circum_ex[0][0] = x_gap;
		_circum_ex[0][1] = y_gap;
		_circum_ex[0][2] = Dist(x_gap, y_gap, 0, _samples[i_rep][0], _samples[i_rep][1], _samples[i_rep][2]);
		_num_circum_ex = 1;

		/*if(min_local_angle<179.5){
		_num_circum_ex=1;
		_circum_ex[0][2]=GetCircumcircle(_samples[i_af],_samples[i_bf],_samples[i_rep],_circum_ex[0][0],_circum_ex[0][1]);
		}else{
		_num_circum_ex=0;
		}*/


#ifdef debug
		if (draw){
			plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, _num_circum_ex, _circum_ex, 0, NULL, 0, NULL, 0);
			plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
		}
#endif
		//the exclusion should maintain the delaunay edges as if we are relocating i_rep

		_num_circum_in = 0;
		_num_line_in = 0;
		n2 = _del[i_rep][_del[i_rep][0]];
		for (V = 1; V <= _del[i_rep][0]; V++){
			n1 = _del[i_rep][V];
			if ((n1 == ip  && n2 == ip1) || (n1 == ip1 && n2 == ip) ||
				(n1 == ip2 && n2 == ip3) || (n1 == ip3 && n2 == ip2)){
				n2 = n1;
				continue;
			}
			ap = FindApex(n1, n2, i_rep, ip, ip1, ip2);

			if (/*!(n1==ip || n2==ip || n1==ip1 || n2==ip1 || n1==ip2 || n2==ip2 || n1==ip3 || n2==ip3)*/
				!(InList(n1, _list) && InList(n2, _list))
				){
				NoThinAngleApex(n1, n2, min_ang, x_gap, y_gap);
				MinAnglePreservMod(n1, n2, min_ang, _line_in[_num_line_in], _samples[i_rep][0], _samples[i_rep][1]);
				_num_line_in++;

				//non-obtuse region 
				_line_in[_num_line_in][0] = _samples[n1][0] - _samples[n2][0];
				_line_in[_num_line_in][1] = _samples[n1][1] - _samples[n2][1];
				_line_in[_num_line_in][2] = -1.0*(_samples[n1][0] * _line_in[_num_line_in][0] + _samples[n1][1] * _line_in[_num_line_in][1]);

				_line_in[_num_line_in][3] = _samples[n2][0] - _samples[n1][0];
				_line_in[_num_line_in][4] = _samples[n2][1] - _samples[n1][1];
				_line_in[_num_line_in][5] = -1.0*(_samples[n2][0] * _line_in[_num_line_in][3] + _samples[n2][1] * _line_in[_num_line_in][4]);
				_num_line_in++;

				//diametral circle 
				_circum_ex[_num_circum_ex][0] = (_samples[n1][0] + _samples[n2][0]) / 2.0;
				_circum_ex[_num_circum_ex][1] = (_samples[n1][1] + _samples[n2][1]) / 2.0;
				_circum_ex[_num_circum_ex][2] = Dist(_circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1], 0, _samples[n1][0], _samples[n1][1], 0);
				_num_circum_ex++;

			}


			if (ap == i_rep){ n2 = n1; continue; }
			_circum_ex[_num_circum_ex][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[ap], _circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1]);
			_num_circum_ex++;
			n2 = n1;
		}

		bool n1_inlist, n2_inlist, ap_inlist;
		n2 = _del[i_rep][_del[i_rep][0]];
		for (V = 1; V <= _del[i_rep][0]; V++){
			ap = _del[i_rep][V];
			if (V == _del[i_rep][0]){ n1 = _del[i_rep][1]; }
			else{ n1 = _del[i_rep][V + 1]; }

			n1_inlist = InList(n1, _list);
			n2_inlist = InList(n2, _list);
			ap_inlist = InList(ap, _list);

			if ((n1_inlist&&n2_inlist) || (ap_inlist&&n2_inlist) || (n1_inlist&&ap_inlist)){ n2 = ap; continue; }
			_circum_in[_num_circum_in][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[ap], _circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1]);
			if (_circum_in[_num_circum_in][2]>2){
				//if it's too large circumcircle, then there is no need to check on it
				//since it already contains the whole domain
				n2 = ap;
				continue;
			}
			//must check if such a circle is empty
			//it could contain other existing sample point so as to it'll already never be formed
			//we check against 1-level delaunay neighbour
			if (!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[n1], n2, ap, i_rep, i_rep, i_rep, i_rep) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[n2], n1, ap, i_rep, i_rep, i_rep, i_rep) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[ap], n2, n1, i_rep, i_rep, i_rep, i_rep) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _list, n1, n2, ap, i_rep, i_rep, i_rep)){
				n2 = ap; continue;
			}
			_num_circum_in++;
			n2 = ap;
		}



#ifdef debug
		if (draw){
			plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, _num_circum_in, _circum_in, 0, NULL, 0, NULL, 0);
			plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
		}
#endif


#ifdef debug
		if (draw){
			plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, _num_circum_ex, _circum_ex, 0, NULL, 0, NULL, 0);
			plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
			n2 = _del[i_rep][_del[i_rep][0]];
			for (V = 1; V <= _del[i_rep][0]; V++){
				n1 = _del[i_rep][V];
				if ((n1 == ip &&n2 == ip1) || (n1 == ip1&&n2 == ip) ||
					(n1 == ip2&&n2 == ip3) || (n1 == ip3&&n2 == ip2)){
					n2 = n1;
					continue;
				}
				ap = FindApex(n1, n2, i_rep, ip, ip1, ip2);
				if (ap == i_rep){ n2 = n1; continue; }
				plot_single_dot_n_layers(layer_plot, ip, _samples, _samples[ap][0], _samples[ap][1], r_plot);
				n2 = n1;
			}
		}
#endif

		//the exclusion region should also maintain the edges of the gap
		//the edges at are not connected to ip
		n2 = _list[_list[0]];
		for (V = 1; V <= _list[0]; V++){
			n1 = _list[V];
			if (n1 == i_rep || n2 == i_rep){
				n2 = n1;
				continue;
			}

			if ((n1 == ip &&n2 == ip1) || (n1 == ip1&&n2 == ip) ||
				(n1 == ip2&&n2 == ip3) || (n1 == ip3&&n2 == ip2)){
				n2 = n1;
				continue;
			}


			if (/*!(n1==ip || n2==ip || n1==ip1 || n2==ip1 || n1==ip2 || n2==ip2 || n1==ip3 || n2==ip3)*/
				!(InList(n1, _list) && InList(n2, _list))
				){
				//diametral circle 
				_circum_ex[_num_circum_ex][0] = (_samples[n1][0] + _samples[n2][0]) / 2.0;
				_circum_ex[_num_circum_ex][1] = (_samples[n1][1] + _samples[n2][1]) / 2.0;
				_circum_ex[_num_circum_ex][2] = Dist(_circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1], 0, _samples[n1][0], _samples[n1][1], 0);
				_num_circum_ex++;
			}

			ap = FindApex(n1, n2, ip, ip1, ip2, ip3);
			if (ap == ip){ n2 = n1; continue; }
			_circum_ex[_num_circum_ex][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[ap], _circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1]);
			_num_circum_ex++;

			n2 = n1;
		}

#ifdef debug
		if (draw){
			plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, _num_circum_ex, _circum_ex, 0, NULL, 0, NULL, 0);
			plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
			for (V = 1; V <= _del[i_rep][0]; V++){
				n1 = _del[i_rep][V];
				ap = FindApex(n1, n2, i_rep, ip, ip1, ip2);
				if (ap == i_rep){ n2 = n1; continue; }
				plot_single_dot_n_layers(layer_plot, ip, _samples, _samples[ap][0], _samples[ap][1], r_plot);
				n2 = n1;
			}
			n2 = _list[_list[0]];
			for (V = 1; V <= _list[0]; V++){
				n1 = _list[V];
				if (n1 == i_rep || n2 == i_rep){
					n2 = n1;
					continue;
				}
				if ((n1 == ip &&n2 == ip1) || (n1 == ip1&&n2 == ip) ||
					(n1 == ip2&&n2 == ip3) || (n1 == ip3&&n2 == ip2)){
					n2 = n1;
					continue;
				}
				ap = FindApex(n1, n2, ip, ip1, ip2, ip3);
				if (ap == ip){ n2 = n1; continue; }
				plot_single_dot_n_layers(layer_plot, ip, _samples, _samples[ap][0], _samples[ap][1], r_plot);
				n2 = n1;

			}
		}
#endif		
		//bulding the grid
		size_t num_active_cells = 16;

		// 4X4 initial cells
		_active_cells_i[0] = 0; _active_cells_j[0] = 0;
		_active_cells_i[1] = 0; _active_cells_j[1] = 1;
		_active_cells_i[2] = 0; _active_cells_j[2] = 2;
		_active_cells_i[3] = 0; _active_cells_j[3] = 3;

		_active_cells_i[4] = 1; _active_cells_j[4] = 0;
		_active_cells_i[5] = 1; _active_cells_j[5] = 1;
		_active_cells_i[6] = 1; _active_cells_j[6] = 2;
		_active_cells_i[7] = 1; _active_cells_j[7] = 3;

		_active_cells_i[8] = 2; _active_cells_j[8] = 0;
		_active_cells_i[9] = 2; _active_cells_j[9] = 1;
		_active_cells_i[10] = 2; _active_cells_j[10] = 2;
		_active_cells_i[11] = 2; _active_cells_j[11] = 3;

		_active_cells_i[12] = 3; _active_cells_j[12] = 0;
		_active_cells_i[13] = 3; _active_cells_j[13] = 1;
		_active_cells_i[14] = 3; _active_cells_j[14] = 2;
		_active_cells_i[15] = 3; _active_cells_j[15] = 3;

		double x_min_grid(10E4), y_min_grid(10E4), x_max_grid(-10E4), y_max_grid(-10E4), s, lx, ly;
		size_t idart, iactive, ii, jj, k;

		for (k = 1; k <= _del[i_rep][0]; k++){
			if (_samples[_del[i_rep][k]][0]<x_min_grid){ x_min_grid = _samples[_del[i_rep][k]][0]; }
			if (_samples[_del[i_rep][k]][1]<y_min_grid){ y_min_grid = _samples[_del[i_rep][k]][1]; }

			if (_samples[_del[i_rep][k]][0]>x_max_grid){ x_max_grid = _samples[_del[i_rep][k]][0]; }
			if (_samples[_del[i_rep][k]][1]>y_max_grid){ y_max_grid = _samples[_del[i_rep][k]][1]; }
		}

		for (k = 1; k <= _list[0]; k++){
			if (_samples[_list[k]][0]<x_min_grid){ x_min_grid = _samples[_list[k]][0]; }
			if (_samples[_list[k]][1]<y_min_grid){ y_min_grid = _samples[_list[k]][1]; }

			if (_samples[_list[k]][0]>x_max_grid){ x_max_grid = _samples[_list[k]][0]; }
			if (_samples[_list[k]][1]>y_max_grid){ y_max_grid = _samples[_list[k]][1]; }
		}

		lx = x_max_grid - x_min_grid;
		ly = y_max_grid - y_min_grid;
		if (lx>ly){ s = lx / 4.0; }
		else{ s = ly / 4.0; }
		double x_new, y_new;


		find = false;
		num_new_samples = 0;

		while (num_new_samples<10){

			RefineGrid(x_min_grid, y_min_grid, s, num_active_cells);
			if (_max_active_cell == num_active_cells){
				break; //unable to find new replacement cuz it's probably lies on a line
			}
#ifdef debug
			if (draw){
				plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, _num_circum_ex, _circum_ex, 0, NULL, 0, NULL, 0);
				plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, _num_circum_in, _circum_in, 0, NULL, 0, NULL, 0);
				plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
				plot_grid_part(layer_plot, ip, _samples, x_min_grid, y_min_grid, s, _active_cells_i, _active_cells_j, num_active_cells, r_plot);
			}
#endif
			if (num_active_cells == 0){
				break;
			}

			for (idart = 0; idart<num_active_cells; idart++){
				iactive = size_t((num_active_cells - 1)*RandNumGenerator());
				ii = _active_cells_i[iactive];
				jj = _active_cells_j[iactive];

				x_new = x_min_grid + (double(ii) + RandNumGenerator())*(s);// getting the position of the dart 
				y_new = y_min_grid + (double(jj) + RandNumGenerator())*(s);

				if (CheckNewReplacement(x_new, y_new)){
					num_new_samples++;
					_samples[i_rep][0] = x_new;
					_samples[i_rep][1] = y_new;

					if (i_rep == ip || i_rep == ip1 || i_rep == ip2 || i_rep == ip3){
						angle = angle_three_points_stathead(_samples[i_bf], _samples[i_rep], _samples[i_intr]);
						angle += angle_three_points_stathead(_samples[i_intr], _samples[i_rep], _samples[i_af]);
					}
					else{
						angle = angle_three_points_stathead(_samples[i_bf], _samples[i_rep], _samples[i_af]);
					}
					if (!find){
						find = true;
						angle_best = angle;
						x_best = x_new;
						y_best = y_new;
					}
					else{
						if (angle<angle_best){
							angle_best = angle;
							x_best = x_new;
							y_best = y_new;
						}
					}

					adjusted = true;
				}
			}
		}

		if (find){
			_samples[i_rep][0] = x_best;
			_samples[i_rep][1] = y_best;
#ifdef debug
			if (draw){
				plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, 0, _circum_in, 0, NULL, 0, NULL, 0);
				plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
			}
#endif	
		}
	}
	return adjusted;
}
bool NonObtuseInjectRepeller(size_t ip, size_t ip1, size_t ip2, bool draw, double min_ang)
{
	//we seek to find a sample to inject so as to destroy the obtuse angle at ip
	//in doing that we remove the edge ip1,ip2, but this we create a new sample that has 4-valent 
	//which is too conservative for non-obtuse
	//thus we remove another edge
	//the candidate for that are ip-ip1, ip-ip2, ip1-ip3 or ip2-ip3	
	//the candidate is discareded if the gap it forms when removed is non-convex
	//Additionally, we try to repel the surrounding samples away to enlarge the resampling space
	//we check all the candidate till we find a successful new replacement 

	double r_plot(0.0080);
	size_t layer_plot(10);

	size_t ip3, n1, V, n2, ap, k, d, bd_flag;
	size_t i;
	bool no_replace, adjusted;
	std::vector <size_t> candid;
	std::vector <size_t> candid_edge;

	double x_mid, y_mid;
	ip3 = FindApex(ip1, ip2, ip, ip, ip, ip);
	if (ip3 == ip){

		//cout<<"Error (0) at NonObtuseInjectRepeller()"<<endl;
		return false;
		
	}


	//candidate ip-ip1
	n1 = FindApex(ip, ip1, ip2, ip2, ip2, ip2);
	if (n1 != ip2){
		candid.push_back(n1);
		candid_edge.push_back(ip);
		candid_edge.push_back(ip1);

		/*angle1=angle_three_points_stathead(_samples[n1],_samples[ip ],_samples[ip2]);
		angle2=angle_three_points_stathead(_samples[n1],_samples[ip1],_samples[ip3]);
		if(angle1<180.0-_tol && angle2<180.0-_tol){
		candid.push_back(n1);
		candid_edge.push_back(ip);
		candid_edge.push_back(ip1);
		}*/
	}
	else{
		//it's a boundary edge
		candid.push_back(ip);
		candid_edge.push_back(ip);
		candid_edge.push_back(ip1);
	}

	//candidate ip-ip2
	n1 = FindApex(ip, ip2, ip1, ip1, ip1, ip1);
	if (n1 != ip1){
		candid.push_back(n1);
		candid_edge.push_back(ip);
		candid_edge.push_back(ip2);

		/*angle1=angle_three_points_stathead(_samples[n1],_samples[ip ],_samples[ip1]);
		angle2=angle_three_points_stathead(_samples[n1],_samples[ip2],_samples[ip3]);
		if(angle1<180.0-_tol && angle2<180.0-_tol){
		candid.push_back(n1);
		candid_edge.push_back(ip);
		candid_edge.push_back(ip2);
		}*/
	}
	else{
		//it's a boundary edge
		candid.push_back(ip);
		candid_edge.push_back(ip);
		candid_edge.push_back(ip2);
	}

	//candidate ip1-ip3
	n1 = FindApex(ip1, ip3, ip2, ip2, ip2, ip2);
	if (n1 != ip2){
		candid.push_back(n1);
		candid_edge.push_back(ip1);
		candid_edge.push_back(ip3);

		/*angle1=angle_three_points_stathead(_samples[n1],_samples[ip1],_samples[ip ]);
		angle2=angle_three_points_stathead(_samples[n1],_samples[ip3],_samples[ip2]);
		if(angle1<180.0-_tol && angle2<180.0-_tol){
		candid.push_back(n1);
		candid_edge.push_back(ip1);
		candid_edge.push_back(ip3);
		}*/
	}
	else{
		//it's a boundary edge
		candid.push_back(ip1);
		candid_edge.push_back(ip1);
		candid_edge.push_back(ip3);
	}

	//candidate ip2-ip3
	n1 = FindApex(ip2, ip3, ip1, ip1, ip1, ip1);
	if (n1 != ip1){
		candid.push_back(n1);
		candid_edge.push_back(ip2);
		candid_edge.push_back(ip3);

		/*angle1=angle_three_points_stathead(_samples[n1],_samples[ip2],_samples[ip ]);
		angle2=angle_three_points_stathead(_samples[n1],_samples[ip3],_samples[ip1]);
		if(angle1<180.0-_tol && angle2<180.0-_tol){
		candid.push_back(n1);
		candid_edge.push_back(ip2);
		candid_edge.push_back(ip3);
		}*/
	}
	else{
		//it's a boundary edge
		candid.push_back(ip2);
		candid_edge.push_back(ip2);
		candid_edge.push_back(ip3);
	}


	//form the available sampling region with each candidate
	//and try to sample


	x_mid = (_samples[ip][0] + _samples[ip1][0] + _samples[ip2][0]) / 3.0;
	y_mid = (_samples[ip][1] + _samples[ip1][1] + _samples[ip2][1]) / 3.0;
	adjusted = false;
	for (V = 0; V<candid.size(); V++){

		if (adjusted){
			//this means that in previous iteration, the samples
			//have been adjusted in NeighbourSamplesRepeller
			//we bring them back to their place
			for (i = 1; i <= _list[0]; i++){
				_samples[_list[i]][0] = _xy_org[i][0];
				_samples[_list[i]][1] = _xy_org[i][1];
			}
		}
		adjusted = false;
		_list[1] = ip;
		_list[2] = ip1;
		_list[3] = ip2;
		_list[4] = ip3;
		if (candid[V] == ip || candid[V] == ip1 || candid[V] == ip2 || candid[V] == ip3){
			_list[0] = 4;
		}
		else{
			_list[0] = 5;
			_list[5] = candid[V];
		}

		SortList(x_mid, y_mid, _list);

#ifdef debug
		if (draw){
			plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, 0, _circum_ex, 0, NULL, 0, NULL, 1);
			plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
		}
#endif

		for (i = 1; i <= _list[0]; i++){
			_xy_org[i][0] = _samples[_list[i]][0];
			_xy_org[i][1] = _samples[_list[i]][1];
		}
		if (NeighbourSamplesRepeller(ip1, ip2, candid_edge[2 * V], candid_edge[2 * V + 1], draw, min_ang, layer_plot, r_plot)){
			adjusted = true;
		}


		//find the neighbouring circumcircles (inclusion regions)
		_num_circum_in = 0;
		bd_flag = 0;
		if (_samples[ip1][2] != 0 && _samples[ip2][2] != 0){
			//ip1-ip2 will be destroyed if new replacment is found
			//restricting the replacement to lie within the region that
			//keep the smooth representation of a boundary edge
			BoundaryEdgesCircles(ip1, ip2, ip1, ip2, ((_samples[ip][0] + _samples[ip1][0]) / 2.0), ((_samples[ip][1] + _samples[ip1][1]) / 2.0));
			bd_flag = size_t(_samples[ip1][2]);
		}
		if (_samples[candid_edge[2 * V]][2] != 0 && _samples[candid_edge[2 * V + 1]][2] != 0 && _list[0] != 4){
			//same thing goes for the candidate edge 
			BoundaryEdgesCircles(candid_edge[2 * V], candid_edge[2 * V + 1], candid_edge[2 * V], candid_edge[2 * V + 1], ((_samples[candid_edge[2 * V]][0] + _samples[candid_edge[2 * V + 1]][0]) / 2.0), ((_samples[candid_edge[2 * V]][1] + _samples[candid_edge[2 * V + 1]][1]) / 2.0));
			bd_flag = size_t(_samples[candid_edge[2 * V]][2]);
		}


		n2 = _list[_list[0]];
		for (i = 1; i <= _list[0]; i++){
			ap = _list[i];
			if (i == _list[0]){ n1 = _list[1]; }
			else{ n1 = _list[i + 1]; }

			NoThinAngleApex(ap, n2, min_ang, x_mid, y_mid);

			_circum_in[_num_circum_in][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[ap], _circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1]);
			if (_circum_in[_num_circum_in][2]>2){
				//if it's too large circumcircle, then there is no need to check on it
				//since it already contains the whole domain
				n2 = ap;
				continue;
			}

			if (!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[n1], n2, ap, _num_samples, _num_samples, _num_samples, _num_samples) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[n2], n1, ap, _num_samples, _num_samples, _num_samples, _num_samples) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[ap], n2, n1, _num_samples, _num_samples, _num_samples, _num_samples) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _list, n1, n2, ap, _num_samples, _num_samples, _num_samples)){
				n2 = ap; continue;
			}
			_num_circum_in++;
			n2 = ap;
		}

#ifdef debug
		if (draw){
			plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, _num_circum_in, _circum_in, 0, NULL, 0, NULL, 1);
			plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
		}
#endif

		no_replace = false;
		if (_num_circum_in>1){
			for (k = 0; k<_num_circum_in - 1; k++){
				for (d = k + 1; d<_num_circum_in; d++){
					if (Dist(_circum_in[k][0], _circum_in[k][1], 0, _circum_in[d][0], _circum_in[d][1], 0)*(1.0 + _tol)*(1.0 + _tol)>_circum_in[k][2] + _circum_in[d][2] + 2.0*sqrt(_circum_in[k][2])*sqrt(_circum_in[d][2])){
						no_replace = true;
						break;
					}
				}
				if (no_replace){ break; }
			}
		}
		if (no_replace){
			continue;
		}

		//find the neighbouring circumcircles (exclusion regions)
		//and find inclusion area to preserve min angle (no thin triangle area)
		_num_circum_ex = 0;
		_num_line_in = 0;
		n2 = _list[_list[0]];
		for (k = 1; k <= _list[0]; k++){
			n1 = _list[k];

			if (((n1 == candid_edge[2 * V] && n2 == candid_edge[2 * V + 1]) || (n2 == candid_edge[2 * V] && n1 == candid_edge[2 * V + 1])) && _list[0] != 4){
				n2 = n1;
				continue;
			}
			ap = FindApex(n1, n2, ip, ip1, ip2, ip3);


			MinAnglePreservMod(n1, n2, min_ang, _line_in[_num_line_in], x_mid, y_mid);
			//plot_single_dot_n_layers(3,ip,_samples,_line_in[_num_line_in][4],_line_in[_num_line_in][5],0.1);
			_num_line_in++;

			//diametral circle 
			_circum_ex[_num_circum_ex][0] = (_samples[n1][0] + _samples[n2][0]) / 2.0;
			_circum_ex[_num_circum_ex][1] = (_samples[n1][1] + _samples[n2][1]) / 2.0;
			_circum_ex[_num_circum_ex][2] = Dist(_circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1], 0, _samples[n1][0], _samples[n1][1], 0);
			_num_circum_ex++;

			//non-obtuse region 
			_line_in[_num_line_in][0] = _samples[n1][0] - _samples[n2][0];
			_line_in[_num_line_in][1] = _samples[n1][1] - _samples[n2][1];
			_line_in[_num_line_in][2] = -1.0*(_samples[n1][0] * _line_in[_num_line_in][0] + _samples[n1][1] * _line_in[_num_line_in][1]);

			_line_in[_num_line_in][3] = _samples[n2][0] - _samples[n1][0];
			_line_in[_num_line_in][4] = _samples[n2][1] - _samples[n1][1];
			_line_in[_num_line_in][5] = -1.0*(_samples[n2][0] * _line_in[_num_line_in][3] + _samples[n2][1] * _line_in[_num_line_in][4]);
			_num_line_in++;

			//neighbour circumcircle
			if (ap == ip || ap == ip1){ n2 = n1; continue; }//in case there is no apex to report		
			_circum_ex[_num_circum_ex][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[ap], _circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1]);
			_num_circum_ex++;
			n2 = n1;
		}

#ifdef debug
		if (draw){
			plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, _num_circum_ex, _circum_ex, 0, NULL, 0, NULL, 1);
			plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
			for (k = 1; k <= _list[0]; k++){
				n1 = _list[k];
				ap = FindApex(n1, n2, ip, ip1, ip2, ip3);
				if (ap == ip || ap == ip1){ n2 = n1; continue; };//in case there is no apex to report
				plot_single_dot_n_layers(layer_plot, ip, _samples, _samples[ap][0], _samples[ap][1], r_plot);
				n2 = n1;
			}
		}
#endif

		//bulding the grid
		size_t num_active_cells = 16;

		// 4X4 initial cells
		_active_cells_i[0] = 0; _active_cells_j[0] = 0;
		_active_cells_i[1] = 0; _active_cells_j[1] = 1;
		_active_cells_i[2] = 0; _active_cells_j[2] = 2;
		_active_cells_i[3] = 0; _active_cells_j[3] = 3;

		_active_cells_i[4] = 1; _active_cells_j[4] = 0;
		_active_cells_i[5] = 1; _active_cells_j[5] = 1;
		_active_cells_i[6] = 1; _active_cells_j[6] = 2;
		_active_cells_i[7] = 1; _active_cells_j[7] = 3;

		_active_cells_i[8] = 2; _active_cells_j[8] = 0;
		_active_cells_i[9] = 2; _active_cells_j[9] = 1;
		_active_cells_i[10] = 2; _active_cells_j[10] = 2;
		_active_cells_i[11] = 2; _active_cells_j[11] = 3;

		_active_cells_i[12] = 3; _active_cells_j[12] = 0;
		_active_cells_i[13] = 3; _active_cells_j[13] = 1;
		_active_cells_i[14] = 3; _active_cells_j[14] = 2;
		_active_cells_i[15] = 3; _active_cells_j[15] = 3;

		double x_min_grid(10E4), y_min_grid(10E4), x_max_grid(-10E4), y_max_grid(-10E4), s, lx, ly;
		size_t idart, iactive, ii, jj;

		for (k = 1; k <= _list[0]; k++){
			if (_samples[_list[k]][0]<x_min_grid){ x_min_grid = _samples[_list[k]][0]; }
			if (_samples[_list[k]][1]<y_min_grid){ y_min_grid = _samples[_list[k]][1]; }

			if (_samples[_list[k]][0]>x_max_grid){ x_max_grid = _samples[_list[k]][0]; }
			if (_samples[_list[k]][1]>y_max_grid){ y_max_grid = _samples[_list[k]][1]; }
		}

		lx = x_max_grid - x_min_grid;
		ly = y_max_grid - y_min_grid;
		if (lx>ly){ s = lx / 4.0; }
		else{ s = ly / 4.0; }
		double x_new, y_new;


		while (true){

#ifdef debug
			if (draw){
				plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, _num_circum_ex, _circum_ex, 0, NULL, 0, NULL, 1);
				plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, _num_circum_in, _circum_in, 0, NULL, 0, NULL, 1);
				plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
				plot_grid_part(layer_plot, ip, _samples, x_min_grid, y_min_grid, s, _active_cells_i, _active_cells_j, num_active_cells, r_plot);
			}
#endif

			RefineGrid(x_min_grid, y_min_grid, s, num_active_cells);
			if (_max_active_cell == num_active_cells){
				break; //unable to find new replacement cuz it's probably lies on a line
			}
#ifdef debug
			if (draw){
				plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, _num_circum_ex, _circum_ex, 0, NULL, 0, NULL, 1);
				plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, _num_circum_in, _circum_in, 0, NULL, 0, NULL, 1);
				plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
				plot_grid_part(layer_plot, ip, _samples, x_min_grid, y_min_grid, s, _active_cells_i, _active_cells_j, num_active_cells, r_plot);
			}
#endif

			if (num_active_cells == 0){
				break;
			}

			for (idart = 0; idart<num_active_cells; idart++){
				iactive = size_t((num_active_cells - 1)*RandNumGenerator());
				ii = _active_cells_i[iactive];
				jj = _active_cells_j[iactive];

				x_new = x_min_grid + (double(ii) + RandNumGenerator())*(s);// getting the position of the dart 
				y_new = y_min_grid + (double(jj) + RandNumGenerator())*(s);


				if (CheckNewReplacement(x_new, y_new)){

					RemoveNodeFromList(_del[ip1], ip2);
					RemoveNodeFromList(_del[ip2], ip1);

					if (_list[0] != 4){
						RemoveNodeFromList(_del[candid_edge[2 * V]], candid_edge[2 * V + 1]);
						RemoveNodeFromList(_del[candid_edge[2 * V + 1]], candid_edge[2 * V]);
					}

					//update delaunay triangulation
					_samples[_num_samples][0] = x_new;
					_samples[_num_samples][1] = y_new;
					_samples[_num_samples][2] = bd_flag;

					for (k = 0; k <= _list[0]; k++){
						_del[_num_samples][k] = _list[k];
					}

					ap = _list[_list[0]];
					for (k = 1; k <= _list[0]; k++){
						if (k == _list[0]){ n1 = _list[1]; }
						else{ n1 = _list[k + 1]; }
						AddNewNodeToDel(_num_samples, _list[k], ap, n1);
						ap = _list[k];
					}
					_num_samples++;
#ifdef debug
					if (draw){
						plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, 0, _circum_in, 0, NULL, 0, NULL, 1);
						plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
						//plot_unit_box(_samples,_del,_num_samples,0.1,0,0,0,1,0,0,0,0,NULL,0,0,NULL);
						//mark_point(_samples[ip][0],_samples[ip][1],_samples,NULL,_num_samples,0,0.03);
					}
#endif					
					return true;

				}
			}
		}
	}

	if (adjusted){
		//this means that in last unsuccessful iteration, the samples
		//have been adjusted in NeighbourSamplesRepeller
		//we bring them back to their place
		for (i = 1; i <= _list[0]; i++){
			_samples[_list[i]][0] = _xy_org[i][0];
			_samples[_list[i]][1] = _xy_org[i][1];
		}
	}
#ifdef debug
	if (draw){
		plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, 0, _circum_ex, 0, NULL, 0, NULL, 1);
		plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
	}
#endif


	return false;
}
bool NonObtuseRelocate(size_t ip, bool draw, double min_ang)
{
	if (_samples[ip][3] != 0){ return false; }

	double r_plot(0.04);
	size_t layer_plot(15);

	size_t V, d, n1, n2, ap, thrd, bd1, bd2, tmp;
	double angle1, angle2, angle3;
	size_t i;

	bool no_replace;
	for (V = 0; V <= _del[ip][0]; V++){
		_list[V] = _del[ip][V];
	}
	_num_circum_in = 0;
	size_t bd_flag = size_t(_samples[ip][2]);
	_num_line_in = 0;
	if (bd_flag != 0){
		//if it's a boundary point 
		if (IsCorner(ip, bd1, bd2)){
			//after sifting few point on the boundary a sample point that is not a corner 
			//becomes corner
			_samples[ip][3] = 1;
			return false;
		}
		//adding the two circles that will maintain the boundary edge 
		//from being distorded when ip moves
		BoundaryEdgesCircles(bd1, bd2, bd1, bd2, _samples[ip][0], _samples[ip][1]);
		/*//put bd1 and bd2 as the start and end of neighbours samples list (_list)
		while( !((_list[1]==bd1 && _list[_list[0]]==bd2)  || (_list[1]==bd2 && _list[_list[0]]==bd1)) ){
		n1=_list[_list[0]];
		for(i=1;i<=_list[0];i++){
		tmp=_list[i];
		_list[i]=n1;
		n1=tmp;
		}
		}*/

		//put bd1 and bd2 as the start and end of neighbours samples list (_list)

		size_t counter(0);
		while ((_samples[_list[1]][2] != bd_flag || _samples[_list[_list[0]]][2] != bd_flag) && (counter<_list[0])) {
			counter++;
			n1 = _list[_list[0]];
			for (i = 1; i <= _list[0]; i++) {
				tmp = _list[i];
				_list[i] = n1;
				n1 = tmp;
			}
		}
		if (_samples[_list[1]][2] != bd_flag || _samples[_list[_list[0]]][2] != bd_flag) {
			//giving up
			std::cout << "Error (0) NonObtuseRelocate2. I AM GIVING UP" << std::endl;
			
			return false;
		}
	}

#ifdef debug
	if (draw){
		plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, 0, NULL, 0, NULL, 0, NULL, 1);
		plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
	}
#endif

	//find the neighbouring circumcircles (inclusion regions)	
	if (_samples[ip][2] == 0){
		n2 = _list[_list[0]];
		for (i = 1; i <= _list[0]; i++){
			ap = _list[i];
			if (i == _list[0]){ n1 = _list[1]; }
			else{ n1 = _list[i + 1]; }

			//first get the circumcircle that prevent thin angle at the apex (new replacement)
			NoThinAngleApex(ap, n2, min_ang, _samples[ip][0], _samples[ip][1]);

			_circum_in[_num_circum_in][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[ap], _circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1]);
			if (_circum_in[_num_circum_in][2]>2){
				//if it's too large circumcircle, then there is no need to check on it
				//since it already contains the whole domain
				n2 = ap;
				continue;
			}
			//must check if such a circle is empty
			//it could contain other existing sample point so as to it'll already never be formed
			//we check against 1-level delaunay neighbour
			if (!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[n1], n2, ap, ip, ip, ip, ip) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[n2], n1, ap, ip, ip, ip, ip) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[ap], n2, n1, ip, ip, ip, ip) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _list, n1, n2, ap, ip, ip, ip)){
				n2 = ap; continue;
			}
			_num_circum_in++;
			n2 = ap;
		}
	}
	else{

		for (i = 1; i< _list[0] - 1; i++){
			n1 = _list[i];
			ap = _list[i + 1];
			n2 = _list[i + 2];

			//first get the circumcircle that prevent thin angle at the apex (new replacement)
			NoThinAngleApex(n1, ap, min_ang, _samples[ip][0], _samples[ip][1]);

			_circum_in[_num_circum_in][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[ap], _circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1]);
			if (_circum_in[_num_circum_in][2]>2){
				//if it's too large circumcircle, then there is no need to check on it
				//since it already contains the whole domain			
				continue;
			}
			//must check if such a circle is empty
			//it could contain other existing sample point so as to it'll already never be formed
			//we check against 1-level delaunay neighbour
			if (!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[n1], n2, ap, ip, ip, ip, ip) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[n2], n1, ap, ip, ip, ip, ip) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _del[ap], n2, n1, ip, ip, ip, ip) ||
				!EmptyCircle(_circum_in[_num_circum_in][0], _circum_in[_num_circum_in][1], _circum_in[_num_circum_in][2], _list, n1, n2, ap, ip, ip, ip)){
				continue;
			}
			_num_circum_in++;
		}
		NoThinAngleApex(ap, n2, min_ang, _samples[ip][0], _samples[ip][1]);
	}

#ifdef debug
	if (draw){
		plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, _num_circum_in, _circum_in, 0, NULL, 0, NULL, 1);
		plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
	}
#endif
	//if there is no area of intersection between the inclusion circle,
	//then for sure there is no replacement
	no_replace = false;
	if (_num_circum_in>1){
		for (V = 0; V<_num_circum_in; V++){
			for (d = V + 1; d<_num_circum_in; d++){
				if (Dist(_circum_in[V][0], _circum_in[V][1], 0, _circum_in[d][0], _circum_in[d][1], 0)*(1.0 + _tol)*(1.0 + _tol)>_circum_in[V][2] + _circum_in[d][2] + 2.0*sqrt(_circum_in[V][2])*sqrt(_circum_in[d][2])){
					no_replace = true;
					break;
				}
			}
			if (no_replace){ break; }
		}
	}
	if (no_replace){
		return false;
	}

	//find the neighbouring circumcircles (exclusion regions)
	//and find inclusion area to preserve the non-obtuse triangle
	_num_circum_ex = 0;

	if (_samples[ip][2] == 0){
		n2 = _list[_list[0]];
		for (V = 1; V <= _list[0]; V++){
			n1 = _list[V];
			if (_list[0] == 3){
				if (V == 3){ thrd = _list[1]; }
				else{ thrd = _list[V + 1]; }
				ap = FindApex(n1, n2, ip, ip, thrd, ip);
			}
			else{
				ap = FindApex(n1, n2, ip, ip, ip, ip);
			}

			MinAnglePreservMod(n1, n2, min_ang, _line_in[_num_line_in], _samples[ip][0], _samples[ip][1]);

			_num_line_in++;

			//diametral circle 
			_circum_ex[_num_circum_ex][0] = (_samples[n1][0] + _samples[n2][0]) / 2.0;
			_circum_ex[_num_circum_ex][1] = (_samples[n1][1] + _samples[n2][1]) / 2.0;
			_circum_ex[_num_circum_ex][2] = Dist(_circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1], 0, _samples[n1][0], _samples[n1][1], 0);
			_num_circum_ex++;

			//non-obtuse region 
			_line_in[_num_line_in][0] = _samples[n1][0] - _samples[n2][0];
			_line_in[_num_line_in][1] = _samples[n1][1] - _samples[n2][1];
			_line_in[_num_line_in][2] = -1.0*(_samples[n1][0] * _line_in[_num_line_in][0] + _samples[n1][1] * _line_in[_num_line_in][1]);

			_line_in[_num_line_in][3] = _samples[n2][0] - _samples[n1][0];
			_line_in[_num_line_in][4] = _samples[n2][1] - _samples[n1][1];
			_line_in[_num_line_in][5] = -1.0*(_samples[n2][0] * _line_in[_num_line_in][3] + _samples[n2][1] * _line_in[_num_line_in][4]);
			_num_line_in++;

			//neighbour circumcircle
			if (ap == ip){ n2 = n1; continue; }//in case there is no apex to report		
			_circum_ex[_num_circum_ex][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[ap], _circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1]);
			_num_circum_ex++;
			n2 = n1;
		}
	}
	else{
		//if it's a boundary sample to relocate
		for (V = 1; V<_list[0]; V++){
			n1 = _list[V];
			n2 = _list[V + 1];
			if (_list[0] == 3){
				if (V == 3){ thrd = _list[1]; }
				else{ thrd = _list[V + 1]; }
				ap = FindApex(n1, n2, ip, ip, thrd, ip);
			}
			else{
				ap = FindApex(n1, n2, ip, ip, ip, ip);
			}

			MinAnglePreservMod(n1, n2, min_ang, _line_in[_num_line_in], _samples[ip][0], _samples[ip][1]);
			_num_line_in++;

			//diametral circle 
			_circum_ex[_num_circum_ex][0] = (_samples[n1][0] + _samples[n2][0]) / 2.0;
			_circum_ex[_num_circum_ex][1] = (_samples[n1][1] + _samples[n2][1]) / 2.0;
			_circum_ex[_num_circum_ex][2] = Dist(_circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1], 0, _samples[n1][0], _samples[n1][1], 0);
			_num_circum_ex++;

			//non-obtuse region 
			_line_in[_num_line_in][0] = _samples[n1][0] - _samples[n2][0];
			_line_in[_num_line_in][1] = _samples[n1][1] - _samples[n2][1];
			_line_in[_num_line_in][2] = -1.0*(_samples[n1][0] * _line_in[_num_line_in][0] + _samples[n1][1] * _line_in[_num_line_in][1]);


			_line_in[_num_line_in][3] = _samples[n2][0] - _samples[n1][0];
			_line_in[_num_line_in][4] = _samples[n2][1] - _samples[n1][1];
			_line_in[_num_line_in][5] = -1.0*(_samples[n2][0] * _line_in[_num_line_in][3] + _samples[n2][1] * _line_in[_num_line_in][4]);
			_num_line_in++;

			//neighbour circumcircle
			if (ap == ip){ n2 = n1; continue; }//in case there is no apex to report		
			_circum_ex[_num_circum_ex][2] = GetCircumcircle(_samples[n1], _samples[n2], _samples[ap], _circum_ex[_num_circum_ex][0], _circum_ex[_num_circum_ex][1]);
			_num_circum_ex++;
		}
	}
#ifdef debug
	if (draw){
		plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, _num_circum_ex, _circum_ex, 0, NULL, 0, NULL, 1);
		plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
		for (V = 1; V <= _list[0]; V++){
			n1 = _list[V];
			if (_list[0] == 3){
				if (V == 3){ thrd = _list[1]; }
				else{ thrd = _list[V + 1]; }
				ap = FindApex(n1, n2, ip, ip, thrd, ip);
			}
			else{
				ap = FindApex(n1, n2, ip, ip, ip, ip);
			}
			if (ap == ip){ n2 = n1; continue; };//in case there is no apex to report
			plot_single_dot_n_layers(layer_plot, ip, _samples, _samples[ap][0], _samples[ap][1], r_plot);
			n2 = n1;
		}
	}
#endif


	//bulding the grid
	size_t num_active_cells = 16;

	// 4X4 initial cells
	_active_cells_i[0] = 0; _active_cells_j[0] = 0;
	_active_cells_i[1] = 0; _active_cells_j[1] = 1;
	_active_cells_i[2] = 0; _active_cells_j[2] = 2;
	_active_cells_i[3] = 0; _active_cells_j[3] = 3;

	_active_cells_i[4] = 1; _active_cells_j[4] = 0;
	_active_cells_i[5] = 1; _active_cells_j[5] = 1;
	_active_cells_i[6] = 1; _active_cells_j[6] = 2;
	_active_cells_i[7] = 1; _active_cells_j[7] = 3;

	_active_cells_i[8] = 2; _active_cells_j[8] = 0;
	_active_cells_i[9] = 2; _active_cells_j[9] = 1;
	_active_cells_i[10] = 2; _active_cells_j[10] = 2;
	_active_cells_i[11] = 2; _active_cells_j[11] = 3;

	_active_cells_i[12] = 3; _active_cells_j[12] = 0;
	_active_cells_i[13] = 3; _active_cells_j[13] = 1;
	_active_cells_i[14] = 3; _active_cells_j[14] = 2;
	_active_cells_i[15] = 3; _active_cells_j[15] = 3;

	double x_min_grid(10E4), y_min_grid(10E4), x_max_grid(-10E4), y_max_grid(-10E4), s, lx, ly, x_new, y_new;
	size_t idart, iactive, ii, jj;

	for (V = 1; V <= _list[0]; V++){
		if (_samples[_list[V]][0]<x_min_grid){ x_min_grid = _samples[_list[V]][0]; }
		if (_samples[_list[V]][1]<y_min_grid){ y_min_grid = _samples[_list[V]][1]; }

		if (_samples[_list[V]][0]>x_max_grid){ x_max_grid = _samples[_list[V]][0]; }
		if (_samples[_list[V]][1]>y_max_grid){ y_max_grid = _samples[_list[V]][1]; }
	}

	lx = x_max_grid - x_min_grid;
	ly = y_max_grid - y_min_grid;
	if (lx>ly){ s = lx / 4.0; }
	else{ s = ly / 4.0; }

	//sampling 
	bool find = false;
	double x_best, y_best, perfect_angle;
	if (_samples[ip][2] == 0){
		perfect_angle = 360.0 / double(_list[0]);
	}
	else{
		perfect_angle = 0;
		for (V = 1; V<_list[0]; V++){
			n1 = _list[V];
			n2 = _list[V + 1];
			perfect_angle += angle_three_points_stathead(_samples[n1], _samples[ip], _samples[n2]);
		}
		perfect_angle /= double(_list[0] - 1);
	}
	angle1 = 360.0;//smallest angle with the best x,y
	size_t num_new_samples(0);


	while (num_new_samples<10){
		RefineGrid(x_min_grid, y_min_grid, s, num_active_cells);
		if (_max_active_cell == num_active_cells){
			break; //unable to find new replacement cuz it's probably lies on a line
		}

#ifdef debug
		if (draw){
			plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, _num_circum_ex, _circum_ex, 0, NULL, 0, NULL, 1);
			plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, _num_circum_in, _circum_in, 0, NULL, 0, NULL, 1);
			plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
			plot_grid_part(layer_plot, ip, _samples, x_min_grid, y_min_grid, s, _active_cells_i, _active_cells_j, num_active_cells, r_plot);
		}
#endif
		if (num_active_cells == 0){
			break;
		}

		for (idart = 0; idart<2.0*num_active_cells; idart++){
			iactive = size_t((num_active_cells - 1)*RandNumGenerator());
			ii = _active_cells_i[iactive];
			jj = _active_cells_j[iactive];

			x_new = x_min_grid + (double(ii) + RandNumGenerator())*(s);// getting the position of the dart 
			y_new = y_min_grid + (double(jj) + RandNumGenerator())*(s);


			//plot_single_dot_n_layers(layer_plot, ip, _samples, x_new, y_new, r_plot);

			if (CheckNewReplacement(x_new, y_new)){
				//ratio1=sqrt(MinToMaxEdgeLengRatio(x_new,y_new));

				num_new_samples++;
				_samples[ip][0] = x_new;
				_samples[ip][1] = y_new;
				angle2 = 360;//smallest angle at the head with this new x,y

				if (_samples[ip][2] == 0){
					n2 = _list[_list[0]];
					for (V = 1; V <= _list[0]; V++){
						n1 = _list[V];
						angle3 = angle_three_points_stathead(_samples[n1], _samples[ip], _samples[n2]);
						if (angle3<angle2){
							angle2 = angle3;
						}
						n2 = n1;
					}
				}
				else{
					for (V = 1; V<_list[0]; V++){
						n1 = _list[V];
						n2 = _list[V + 1];
						angle3 = angle_three_points_stathead(_samples[n1], _samples[ip], _samples[n2]);
						if (angle3<angle2){
							angle2 = angle3;
						}
					}
				}



				if (!find){
					//if it's the first one
					find = true;
					x_best = x_new;
					y_best = y_new;
					angle1 = angle2;

				}
				else{
#ifdef debug
					if (_samples[ip][2] == 0){
						if (perfect_angle - angle2 + _smoothness_factor<0.0 || perfect_angle - angle1 + _smoothness_factor<0){
							std::cout << "Error (0) at NonObtuseRelocate2()" << std::endl;
							
						}
					}
#endif

					if (abs(perfect_angle - angle2)<abs(perfect_angle - angle1)){
						//if(factor<best){
						x_best = x_new;
						y_best = y_new;
						angle1 = angle2;
						//best=factor;
					}
				}
			}
		}
	}

	if (find){
		//update delaunay triangulation
		//move ip to the new location, connect it with the neighbour nodes
		//CheckNewReplacement(x_best,y_best);

		_samples[ip][0] = x_best;
		_samples[ip][1] = y_best;
		for (V = 1; V <= _del[ip][0]; V++){
			RemoveNodeFromList(_del[_del[ip][V]], ip);
		}
		for (V = 0; V <= _list[0]; V++){
			_del[ip][V] = _list[V];
		}
		//ap=_list[_list[0]];
		for (V = 1; V <= _list[0]; V++){
			//if(V==_list[0]){n1=_list[1];}
			//else{n1=_list[V+1];}
			AddNewNodeToDel(ip, _list[V], ap, n1);
			//ap=_list[V];
		}

#ifdef debug
		if (draw){
			plot_n_layers(layer_plot, ip, _samples, _del, _num_samples, r_plot, 1, 0, 1, 0, NULL, 0, 0, _circum_in, 0, NULL, 0, NULL, 1);
			plot_dots_n_layers_orange_list(layer_plot, ip, _samples, _list, r_plot);
		}
#endif			
		return true;
	}


	return false;
}

bool NonObtuseRelocate_Caller(size_t ip, double min_ang)
{
	//pick the right vertex to relocate in case of corner vertex

	if (_samples[ip][3] != 0){
		size_t n1, n2;
		bool obt = ObtuseHead(ip, n1, n2);
		//if it's a corner, we can't relocate it cuz this will damage a sharp feature
		//thus we relocate one of the other obtuse triangle heads
		//preferably an internal point, not a point on an edge
		if (_samples[n2][2] == 0){
			//if n2 is interior
			std::swap(ip, n2);
		}
		else if (_samples[n1][2] == 0){
			//if n1 is interior
			std::swap(ip, n1);
		}
		else if (_samples[n1][3] == 0){
			//if n1 is not a corner .i.e an edge point 
			std::swap(ip, n1);
		}
		else if (_samples[n2][3] == 0){
			//if n2 is not a corner .i.e an edge point 
			std::swap(ip, n2);
		}
		else if (obt){
			std::cout << "Error (0) at NonObtuseRelocateAttemp. Do not know what to do!!!" << std::endl;
		}
		else { return false; } //this if all point are coners and it's not an obtuse angle already
	}

	return NonObtuseRelocate(ip, 0, min_ang);
}
bool NonObtuseSift_Caller(size_t ip, size_t n1, size_t n2, double min_ang)
{
	//pick the right vertex to sift 

	if (_samples[ip][3] != 0){
		//if it's a corner, we can't relocate it cuz this will damage a sharp feature
		//thus we relocate one of the other obtuse triangle heads
		//preferably an internal point, not a point on an edge
		if (_samples[n2][2] == 0){
			//if n2 is interior
			std::swap(ip, n2);
		}
		else if (_samples[n1][2] == 0){
			//if n1 is interior
			std::swap(ip, n1);
		}
		else if (_samples[n1][3] == 0){
			//if n1 is not a corner .i.e an edge point 
			std::swap(ip, n1);
		}
		else if (_samples[n2][3] == 0){
			//if n2 is not a corner .i.e an edge point 
			std::swap(ip, n2);
		}
		else{
			std::cout << "Error (1) at NonObtuseReTriangulation2. Do not know what to do!!!" << std::endl;
			
		}
	}
	for (size_t d = 1; d <= _del[ip][0]; d++){
		size_t ip1 = _del[ip][d];
		if (_samples[ip1][3] != 0){ continue; }//if it's a corner			
		if (NonObtuseSift(ip, ip1, 0, min_ang)){
			return true;
		}
	}
	return false;
}
bool NonObtuseSiftAttractor_Caller(size_t ip, size_t n1, size_t n2, double min_ang)
{
	//pick the right vertex to sift 

	if (_samples[ip][3] != 0){
		//if it's a corner, we can't relocate it cuz this will damage a sharp feature
		//thus we relocate one of the other obtuse triangle heads
		//preferably an internal point, not a point on an edge
		if (_samples[n2][2] == 0){
			//if n2 is interior
			std::swap(ip, n2);
		}
		else if (_samples[n1][2] == 0){
			//if n1 is interior
			std::swap(ip, n1);
		}
		else if (_samples[n1][3] == 0){
			//if n1 is not a corner .i.e an edge point 
			std::swap(ip, n1);
		}
		else if (_samples[n2][3] == 0){
			//if n2 is not a corner .i.e an edge point 
			std::swap(ip, n2);
		}
		else{
			std::cout << "Error (2) at NonObtuseReTriangulation2. Do not know what to do!!!" << std::endl;
		}
	}

	for (size_t d = 1; d <= _del[ip][0]; d++){
		size_t ip1 = _del[ip][d];
		if (_samples[ip1][3] != 0){ continue; }//if it's a corner			
		if (NonObtuseSiftAttractor(ip, ip1, 0, min_ang)){
			return true;
		}
	}
	return false;
}

void ResetMesh(char*NodeFileName, char*EleFileName)
{
	//std::cout <<"RESETTING" << std::endl;
	//reset everything by re-reading everything 
	ReadTriangleInput(NodeFileName, EleFileName);

}
void NonObtuseReTriangulation(double min_ang)
{

	
	size_t V, ittr(0), num_obtuse, n1, n2;
	std::vector<size_t> sample_id;
	double ratio;
	ratio = 0;

	

	while (abs(ratio - 100)>_tol){
		ittr++;

		//*************************** Relocation 
		for (V = 0; V<_num_samples; V++){
			NonObtuseRelocate_Caller(V, min_ang);
		}



		//*************************** Sifting
		for (V = 0; V<_num_samples; V++){
			if (!ObtuseHead(V, n1, n2)){ continue; }
			if (NonObtuseSift_Caller(V, n1, n2, min_ang)){
				_numSifting++;
			}
			else{
				_numSifting_failed++;
			}
		}



		//*************************** Sifting Attractor
		for (V = 0; V<_num_samples; V++){
			if (!ObtuseHead(V, n1, n2)){ continue; }
			if (NonObtuseSiftAttractor_Caller(V, n1, n2, min_ang)){
				_numSiftingAtt++;
			}
			else{
				_numSiftingAtt_failed++;
			}
		}



		//*************************** Injection
		for (V = 0; V<_num_samples; V++){
			if (!ObtuseHead(V, n1, n2)){ continue; }
			if (NonObtuseInjection(V, n1, n2, 0, min_ang)){
				_numInject++;
			}
			else{
				_numInject_failed++;
			}
		}

		//*************************** Injection Repeller
		for (V = 0; V<_num_samples; V++){
			if (!ObtuseHead(V, n1, n2)){ continue; }
			if (NonObtuseInjectRepeller(V, n1, n2, 0, min_ang)){
				_numInjectRep++;
			}
			else{
				_numInjectRep_failed++;
			}
		}


#ifdef debug			
		//if (false){
		//	plot_unit_box(_samples, _del, _num_samples, 0.1, 1, 0, 0, 1, 0, 0, 0, 0, NULL, 0, 0, NULL);
		//}
		VerifyMeshNonObtuse(min_ang, num_obtuse);
		ratio = 100.0*(double(_num_obtuse_input) - double(num_obtuse)) / double(_num_obtuse_input);
		printf("ittr = %ld, _num_samples = %ld, Reduction ratio= %lf, num_obtuse=%ld\n", ittr, _num_samples, ratio, num_obtuse);
#endif
	}
	printf("\n Sifting=%ld, Sifting_Att=%ld, Inject=%ld, Inject_Rep=%ld\n", _numSifting, _numSiftingAtt, _numInject, _numInjectRep);
	printf("\n Sifting_failed=%ld, Sifting_Att_failed=%ld, Inject_failed=%ld, Inject_Rep_failed=%ld\n", _numSifting_failed, _numSiftingAtt_failed, _numInject_failed, _numInjectRep_failed);
}


//int main(int argc, char **argv)
//{
int MeshOpt2D(char *NodeFileName, char* EleFileName)
{
	
	//srand (time(NULL) ); // activate for different experiments
	InitiateRandNumGenerator(rand()); //Initialize 	
	

	ReadTriangleInput(NodeFileName, EleFileName);
	VerifyMeshNonObtuse(_min_angle_input, _num_obtuse_input);
	
	

	std::cout << "Input number of samples = " << _num_samples << std::endl;
	std::cout << "Input Min angle = " << _min_angle_input << std::endl;
	std::cout << "Input Max angle = " << _max_angle_input << std::endl;
	std::cout << "Input number of obtuse angles = " << _num_obtuse_input << std::endl;

	std::cout << "\n \n             **********Non-Obtuse Re-Triangulation Starts**********" << std::endl;
	//WriteMesh(argv[1], argv[2], _num_samples, _min_angle_input, _max_angle_input);	
	NonObtuseReTriangulation(_min_angle_input);
	std::cout << "            **********Non-Obtuse Re-Triangulation Ends***********" << std::endl;



	//*******Reporting*******////
	double max_ang, min_ang;
	get_max_min_angles(_num_samples, _del, _samples, max_ang, min_ang, FindApex);
	std::cout << "\n" << std::endl;


	size_t num_obtuse;
	VerifyMeshNonObtuse(_min_angle_input, num_obtuse);

	std::cout << "Output number of samples = " << _num_samples << std::endl;
	std::cout << "Output Min angle = " << min_ang << std::endl;
	std::cout << "Output Max angle = " << max_ang << std::endl;
	std::cout << "Output number of obtuse angles = " << num_obtuse << std::endl;

	//plot_unit_box(_samples, _del, _num_samples, 0.02, 1, 0, 0, 1, 0, 0, 0, 0, NULL, 0, 0, NULL);

	//WriteMesh(argv[1], argv[2], _num_samples, min_ang, max_ang);


	std::cout << "\n\nMission Accomplished" << std::endl;

	return 0;
}