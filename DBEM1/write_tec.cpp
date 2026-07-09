#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "Geo.h"

using namespace std;
static int write_header(int num_result_per_node, int plot_result, ofstream & out);
static int write_zone(int num_total_node, int num_ele, int num_node_per_ele, int num_result_per_node, Point * node_array, int ** node_list, double ** node_result, Point * center_node, double ** center_result, char* zone_title, int plot_result, ofstream & out);
static int cal_center(int num_ele, int num_result_per_node, Point * node_array, int ** node_list, double ** node_result, Point * center_node, double ** center_result);
/*!\brief write a data file for plotting variables on surfaces in tecplot
 */
int write_tecplot(int num_total_node, int num_ele, int num_node_per_ele, int num_result_per_node, Point * node_array, int ** node_list, double ** node_result, int plot_result, int write_disp, char * out1, double amplitude)
{

	int i, j;

	Point * center_node = 0;
	double ** center_result = 0;
	ofstream out(out1);

	Cube m_cube;

	write_header(num_result_per_node, plot_result, out);

	center_node = new Point[num_ele];
	center_result = new double *[num_ele];

	for (i = 0; i < num_ele; ++i)
		center_result[i] = new double[num_result_per_node];

	if (num_node_per_ele == 8) //calculate the center point
	{

		cal_center(num_ele, num_result_per_node, node_array, node_list, node_result, center_node, center_result);

	} //end if num_node_per_ele

	char zone_title_temp1[12] = { "original" };

	write_zone(num_total_node, num_ele, num_node_per_ele, num_result_per_node, node_array, node_list, node_result, center_node, center_result, zone_title_temp1, plot_result, out);


	double amplify = 1.0;

	double * max_array = new double[num_total_node];

	double max_result, max_coord;

	// 确定变形最大值
	max_result = 0.0;
	for (i = 0; i < num_total_node; ++i)
	{
		for (j = 0; j < DIMENSION; ++j)
		{
			if (max_result < fabs(node_result[i][j]))
				max_result = fabs(node_result[i][j]);
		}
	}

	// 确定求解域尺寸
	AssignBox(m_cube, node_array, num_total_node);
	max_coord = m_cube.length;

	//cout << "amplify=" << amplify << endl;


	Point * deformed_node_array;

	Point * deformed_center_node;

	deformed_node_array = new Point[num_total_node];

	deformed_center_node = new Point[num_ele];

	if ((write_disp == 1))
	{


		for (i = 0; i < num_total_node; ++i)
			for (j = 0; j < DIMENSION; ++j)
				deformed_node_array[i].pt[j] = node_array[i].pt[j] + amplify * node_result[i][j];

		if (num_node_per_ele == 8)
		{

			for (i = 0; i < num_ele; ++i)
				for (j = 0; j < DIMENSION; ++j)
					deformed_center_node[i].pt[j] = center_node[i].pt[j] + amplify * center_result[i][j];
		}

		char zone_title_temp2[12] = { "deformed" };
		write_zone(num_total_node, num_ele, num_node_per_ele, num_result_per_node, deformed_node_array, node_list, node_result, deformed_center_node, center_result, zone_title_temp2, plot_result, out);


	}



	out.close();

	if (center_node)
		delete[] center_node;

	if (center_result)
	{

		for (i = 0; i < num_ele; ++i)
			delete[] center_result[i];

		delete[] center_result;
	}


	return 0;

}

static int write_header(int num_result_per_node, int plot_result, ofstream & out)
{
	int i;
	out << "TITLE = \"BEM: 3D SURFACE ZONES\"" << endl;
	out << "VARIABLES = \"X\", \"Y\", \"Z\"";

	if (plot_result == 1)
		for (i = 0; i < num_result_per_node; ++i)
		{
			out << ", \"var" << i << "\"";
		}

	out << endl;
	return 0;

}

static int cal_center(int num_ele, int num_result_per_node, Point * node_array, int ** node_list, double ** node_result, Point * center_node, double ** center_result)
{
	int i, j, k;


	for (i = 0; i < num_ele; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			center_node[i].pt[j] = 0.0;

			for (k = 0; k < 4; ++k)
				center_node[i].pt[j] += -0.25 * node_array[node_list[i][k] - 1].pt[j];

			for (k = 4; k < 8; ++k)
				center_node[i].pt[j] += 0.5 * node_array[node_list[i][k] - 1].pt[j];

		}

		for (j = 0; j < num_result_per_node; ++j)
		{
			center_result[i][j] = 0.0;

			for (k = 0; k < 4; ++k)
				center_result[i][j] += -0.25 * node_result[node_list[i][k] - 1][j];

			for (k = 4; k < 8; ++k)
				center_result[i][j] += 0.5 * node_result[node_list[i][k] - 1][j];
		}


	}

	return 0;

}

static int write_zone(int num_total_node, int num_ele, int num_node_per_ele, int num_result_per_node, Point * node_array, int ** node_list, double ** node_result, Point * center_node, double ** center_result, char* zone_title, int plot_result, ofstream & out)
{
	int i, j;

	switch (num_node_per_ele)
	{

	case 3:
		out << "ZONE " << "T=" << "\"" << zone_title << "\"" << ", N=" << num_total_node << ", E=" << num_ele << ", F=FEPOINT, ET=TRIANGLE" << endl;
		break;

	case 4:
		out << "ZONE " << "T=" << "\"" << zone_title << "\"" << ", N=" << num_total_node << ", E=" << num_ele << ", F=FEPOINT, ET=QUADRILATERAL" << endl;
		break;

	case 6:
		out << "ZONE " << "T=" << "\"" << zone_title << "\"" << ", N=" << num_total_node << ", E=" << 4 * num_ele << ", F=FEPOINT, ET=TRIANGLE" << endl;
		break;

	case 8:
		out << "ZONE " << "T=" << "\"" << zone_title << "\"" << ", N=" << num_total_node + num_ele << ", E=" << 4 * num_ele << ", F=FEPOINT, ET=QUADRILATERAL" << endl;
		break;

	case 9:
		out << "ZONE " << "T=" << "\"" << zone_title << "\"" << ", N=" << num_total_node << ", E=" << 4 * num_ele << ", F=FEPOINT, ET=QUADRILATERAL" << endl;
		break;



	}

	//print node coordinates
	int prec = 10;


	out.precision(prec);

	for (i = 0; i < num_total_node; ++i)
	{
		for (j = 0; j < 3; ++j)
			out << node_array[i].pt[j] << '\t';

		if (plot_result == 1)
			for (j = 0; j < num_result_per_node; ++j)
				out << node_result[i][j] << '\t';

		out << endl;
	}

	if (num_node_per_ele == 8)
	{
		for (i = 0; i < num_ele; ++i)
		{
			for (j = 0; j < 3; ++j)
				out << center_node[i].pt[j] << '\t';

			if (plot_result == 1)
				for (j = 0; j < num_result_per_node; ++j)
					out << center_result[i][j] << '\t';

			out << endl;
		}

	}

	//print node list

	switch (num_node_per_ele)
	{

	case 3:

		for (i = 0; i < num_ele; ++i)
		{
			for (j = 0; j < num_node_per_ele; ++j)
				out << node_list[i][j] << '\t';

			out << endl;
		}

		break;

	case 4:

		for (i = 0; i < num_ele; ++i)
		{
			for (j = 0; j < num_node_per_ele; ++j)
				out << node_list[i][j] << '\t';

			out << endl;
		}

		break;

	case 6:

		for (i = 0; i < num_ele; ++i)
		{
			out << node_list[i][0] << '\t' << node_list[i][3] << '\t' << node_list[i][5] << endl;
			out << node_list[i][3] << '\t' << node_list[i][4] << '\t' << node_list[i][5] << endl;
			out << node_list[i][3] << '\t' << node_list[i][1] << '\t' << node_list[i][4] << endl;
			out << node_list[i][5] << '\t' << node_list[i][4] << '\t' << node_list[i][2] << endl;
		}

		break;

	case 8:

		for (i = 0; i < num_ele; ++i)
		{
			out << node_list[i][0] << '\t' << node_list[i][4] << '\t' << num_total_node + i + 1 << '\t' << node_list[i][7] << endl;
			out << node_list[i][4] << '\t' << node_list[i][1] << '\t' << node_list[i][5] << '\t' << num_total_node + i + 1 << endl;
			out << num_total_node + i + 1 << '\t' << node_list[i][5] << '\t' << node_list[i][2] << '\t' << node_list[i][6] << endl;
			out << node_list[i][7] << '\t' << num_total_node + i + 1 << '\t' << node_list[i][6] << '\t' << node_list[i][3] << endl;
		}

		break;

	case 9:

		for (i = 0; i < num_ele; ++i)
		{
			out << '\t' << node_list[i][0] << '\t' << node_list[i][4] << '\t' << node_list[i][8] << '\t' << node_list[i][7] << endl;
			out << '\t' << node_list[i][4] << '\t' << node_list[i][1] << '\t' << node_list[i][5] << '\t' << node_list[i][8] << endl;
			out << '\t' << node_list[i][8] << '\t' << node_list[i][5] << '\t' << node_list[i][2] << '\t' << node_list[i][6] << endl;
			out << '\t' << node_list[i][7] << '\t' << node_list[i][8] << '\t' << node_list[i][6] << '\t' << node_list[i][3] << endl;

		}

		break;



	}


	return 0;
}



