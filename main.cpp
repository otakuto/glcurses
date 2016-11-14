#include <iostream>
#include <vector>
#include <cmath>
#include <ncurses.h>
#include <unistd.h>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

int width;
int height;

std::vector<std::vector<char>> buffer;

void drawline(int x0, int y0, int x1, int y1)
{
	int dx = std::abs(x1 - x0);
	int dy = std::abs(y1 - y0);
	int sx = x0 < x1 ? 1 : -1;
	int sy = y0 < y1 ? 1 : -1;
	int e = dx - dy;

	auto f = [&](int x, int y, char v)
	{
		if (0 <= x && x < width && 0 <= y && y < height)
		{
			buffer[y][x] = v;
		}
	};

	while (true)
	{
		//buffer[y0][x0] = ' ';
		f(x0,y0,' ');
		if (x0 == x1 && y0 == y1)
		{
			return;
		}

		int e2 = e * 2;
		if (e2 > -dy)
		{
			e -= dy;
			x0 += sx;
		}
		if (e2 < dx)
		{
			e += dx;
			y0 += sy;
		}
	}
}

void drawtri(int x1, int y1, int x2, int y2, int x3, int y3)
{
	drawline(x1, y1, x2, y2);
	drawline(x2, y2, x3, y3);
	drawline(x3, y3, x1, y1);
}

void buffer_init()
{
	buffer.resize(height);
	for (int i = 0; i < height; ++i)
	{
		buffer[i].resize(width);
	}
}

void buffer_clear()
{
	for (int i = 0; i < height; ++i)
	{
		for (int j = 0; j < width; ++j)
		{
			buffer[i][j] = '\0';
		}
	}
}

Eigen::Matrix<double, 4, 4> viewport(int offX, int offY, double width, double height, double n = 0, double f = 1)
{
	Eigen::Matrix<double, 4, 4> matrix;
	matrix <<
		width / 2, 0, 0, offX + (width / 2),
		0, height / 2, 0, offY + (height / 2),
		0, 0, (f - n) / 2, (f + n) / 2,
		0, 0, 0, 1;
	return matrix;
}

Eigen::Matrix<double, 4, 4> perspective(double angle, double hw, double fz, double nz)
{
	Eigen::Matrix<double, 4, 4> matrix;
	matrix <<
		1 / (std::tan(angle / 2) * hw), 0, 0, 0,
		0, 1 / std::tan(angle / 2), 0, 0,
		0, 0, fz*(1/(fz-nz)), 1,
		0, 0, fz*(-nz/(fz-nz)), 0;
	return matrix;
}

Eigen::Matrix<double, 4, 4> lookAt(Eigen::Vector3d eye, Eigen::Vector3d center, Eigen::Vector3d up)
{
	Eigen::Matrix<double, 4, 4> matrix = Eigen::Matrix4d::Identity();
	Eigen::Vector3d forward = center - eye;

	forward.normalize();

	auto side = forward.cross(up);
	side.normalize();

	up = side.cross(forward);

	matrix(0, 0) = side[0];
	matrix(1, 0) = side[1];
	matrix(2, 0) = side[2];

	matrix(0, 1) = up[0];
	matrix(1, 1) = up[1];
	matrix(2, 1) = up[2];

	matrix(0, 2) = -forward[0];
	matrix(1, 2) = -forward[1];
	matrix(2, 2) = -forward[2];

	Eigen::Affine3d t(Eigen::Translation<double, 3>(-eye[0], -eye[1], -eye[2]));
	return matrix * t.matrix();
}

int main()
{
	initscr();
	clear();
	noecho();
	nonl();
	cbreak();
	curs_set(0);
	//halfdelay(1);
	timeout(1);
	nodelay(stdscr, true);
	keypad(stdscr, true);
	getmaxyx(stdscr, height, width);
	notimeout(stdscr, false);

	buffer_init();

	const double pi = std::acos(-1);
	/*
	start_color();
	init_pair(2, COLOR_BLACK, COLOR_WHITE);
	attron(COLOR_PAIR(2));
	*/
	int time = 0;
	auto vp = viewport(0, 0, width-1, height-1);
	auto per = perspective(pi / 2, static_cast<double>(width) / static_cast<double>(height), 1, 1024);
	while (true)
	{
		auto eye = Eigen::Vector4d(2,0,2,0);
		auto look = lookAt({eye[0],eye[1],eye[2]},{0,0,0},{0,1,0});
		//auto look = lookAt({eye[0],eye[1],eye[2]},{std::sin((pi / 1024) * time),0,std::cos((pi / 1024) * time)},{0,1,0});
		usleep(1666);
		clear();
		buffer_clear();
		//auto f = (time % 1000) * 0.01;
		auto f = -2;
		Eigen::Vector3d q1(1,0,f);
		Eigen::Vector3d q2(0,1,f);
		Eigen::Vector3d q3(-1,0,f);
		Eigen::Vector4d p1;
		Eigen::Vector4d p2;
		Eigen::Vector4d p3;
		//auto ro = Eigen::AngleAxisd((pi / 1024) * time, Eigen::Vector3d::UnitY());
		//q1 = ro * q1;
		//q2 = ro * q2;
		//q3 = ro * q3;
		p1 = {q1[0], q1[1], q1[2], 1};
		p2 = {q2[0], q2[1], q2[2], 1};
		p3 = {q3[0], q3[1], q3[2], 1};
		//p1 = vp * per * look * (p1 / f);
		//p2 = vp * per * look * (p2 / f);
		//p3 = vp * per * look * (p3 / f);
		p1 = vp * per * (p1 / f);
		p2 = vp * per * (p2 / f);
		p3 = vp * per * (p3 / f);
		drawtri(p1[0],p1[1],p2[0],p2[1],p3[0],p3[1]);

		if(false){
			Eigen::Vector2d v1(-10,-10);
			Eigen::Vector2d v2(-10,10);
			Eigen::Vector2d v3(10,-10);
			//Eigen::Matrix3d rot;
			auto rot = Eigen::Rotation2D<double>((pi / 1024) * time);
			v1 = rot * v1;
			v2 = rot * v2;
			v3 = rot * v3;
			v1[0] += 20;
			v1[1] += 20;
			v2[0] += 20;
			v2[1] += 20;
			v3[0] += 20;
			v3[1] += 20;
			//std::cerr << rot << std::endl;
			//std::cerr <<  << std::endl;
			//std::cerr << "v1: " << std::endl << v1 << std::endl;
			//std::cerr << "v2: " << std::endl << v2 << std::endl;
			//std::cerr << "v3: " << std::endl << v3 << std::endl;
			drawtri(v1[0],v1[1],v2[0],v2[1],v3[0],v3[1]);
		}
		++time;
		for (int i = 0; i < height; ++i)
		{
			for (int j = 0; j < width; ++j)
			{
				if (buffer[i][j] == ' ')
				{
					mvaddch(i,j,' ' | A_REVERSE);
				}
			}
		}
		refresh();
	}
	endwin();
}

