#include <future> 
#include <chrono>
#include <iomanip>
#include <filesystem>
#include <iostream>
#include <algorithm> 
#include <vector>
#include <string>
#include <limits>
#include <Windows.h>
#include <sstream>

#include "rang.hpp"
#include "OCR_font_STL.h"

#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/bounding_box.h>



#define M_PI 3.14159265358979323846

//<< Red << << ColorEnd <<
auto ColorEnd = [](std::ostream& os) -> std::ostream& { return os << rang::fg::reset; };
auto Red = [](std::ostream& os) -> std::ostream& { return os << rang::fg::red; };
auto Green = [](std::ostream& os) -> std::ostream& { return os << rang::fg::green; };
auto Yellow = [](std::ostream& os) -> std::ostream& { return os << rang::fg::yellow; };
auto Blue = [](std::ostream& os) -> std::ostream& { return os << rang::fg::blue; };
auto Magenta = [](std::ostream& os) -> std::ostream& { return os << rang::fg::magenta; };
auto Cyan = [](std::ostream& os) -> std::ostream& { return os << rang::fg::cyan; };
auto Gray = [](std::ostream& os) -> std::ostream& { return os << rang::fg::gray; };


namespace PMP = CGAL::Polygon_mesh_processing;
namespace fs = std::filesystem;
bool DEBUG = false;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

void resizeConsole(SHORT width, SHORT height) {
	HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
	SMALL_RECT windowSize = { 0, 0, width - 1, height - 1 }; // Minus 1 because it is zero-based

	// Set the console window size
	SetConsoleWindowInfo(hStdout, TRUE, &windowSize);

	// Update the buffer size to at least as large as the window
	COORD bufferSize = { static_cast<SHORT>(width), static_cast<SHORT>(height) };
	SetConsoleScreenBufferSize(hStdout, bufferSize);
}

void get_dimensions(Mesh mesh, double& modelWidth, double& modelLength, double& modelHeight) {
	std::vector<Point> points;
	for (auto v : mesh.vertices()) {
		points.push_back(mesh.point(v));
	}
	Kernel::Iso_cuboid_3 bbox = CGAL::bounding_box(points.begin(), points.end());
	modelWidth = static_cast<double>(bbox.xmax() - bbox.xmin());
	modelLength = static_cast<double>(bbox.ymax() - bbox.ymin());
	modelHeight = static_cast<double>(bbox.zmax() - bbox.zmin());
	if (DEBUG) std::cout << Yellow << "      Dimensions:" << ColorEnd
		<< "  (W"
		<< modelWidth << "  L"
		<< modelLength << "  H"
		<< modelHeight << ")" << std::endl;
}

void scaleMesh(Mesh& mesh, double XYscale, double XYtopscale, double Zscale, double zThreshold) {
	for (auto v : mesh.vertices()) {
		Point& point = mesh.point(v);
		double new_x, new_y, new_z;

		if (point.z() > zThreshold) {
			new_x = point.x() * XYtopscale;
			new_y = point.y() * XYtopscale;
		}
		else {
			new_x = point.x() * XYscale;
			new_y = point.y() * XYscale;
		}
		new_z = point.z() * Zscale;
		mesh.point(v) = Point(new_x, new_y, new_z);
	}
}

void translate_mesh(Mesh& mesh, const Vector& translation_vector) {
	if (DEBUG) std::cout << Yellow << "      Applying translation:  " << ColorEnd << translation_vector << std::endl;
	for (auto v : mesh.vertices()) {
		mesh.point(v) = mesh.point(v) + translation_vector;
	}
}

bool write_STL(const std::string& filename, const Mesh& mesh) {
	fs::path filepath(filename);
	if (DEBUG) std::cout << Yellow << "      Writting STL file:  " << ColorEnd << filepath.filename() << std::endl;
	if (!CGAL::IO::write_polygon_mesh(filename, mesh, CGAL::parameters::stream_precision(10))) {
		std::cerr << Red << "Error: Cannot write the STL file:  " << ColorEnd << filepath.filename() << std::endl;
		return false;
	}
	return true;
}


bool read_STL_data(const std::string& identifier, Mesh& mesh) {
	mesh.clear();
	for (const auto& data : FONT_STL) {
		if (data.key == identifier) { // Convert char to string for comparison
			if (DEBUG) std::cout << Yellow << "      Reading STL Data:  " << ColorEnd << identifier << std::endl;
			std::istringstream iss(std::string(reinterpret_cast<const char*>(data.data), data.size), std::ios::binary);
			if (CGAL::IO::read_STL(iss, mesh)) { // Ensure this matches the actual function available in CGAL
				return true;
			}
			break;
		}
	}
	std::cerr << Red << "      Error: No STL data available for:  " << ColorEnd << identifier << std::endl;
	return false;
}

void create_fixture(std::string ID_Str, Mesh Fixture_Mesh, Mesh& Result_Mesh) {
	bool lastWasDigit = false;
	double offsetX = -6.5, offsetY = -7.5, offsetZ = 4.0;
	double XYscale = 0.18, XYtopscale = 0.18, Zscale = 0.30;
	double zThreshold = 0.1;
	double Xspacing = 0.8, Yspacing = 2.9;
	double zDepth = -1.0;
	Mesh Tag_Mesh;

	std::transform(ID_Str.begin(), ID_Str.end(), ID_Str.begin(), [](unsigned char c) { return std::toupper(c); });

	for (char c : ID_Str) {
		Mesh Letter_Mesh;
		double FontWidth = 0.0, FontLength = 0.0, FontHeight = 0.0;

		if (!read_STL_data(std::string(1, c), Letter_Mesh)) continue;

		get_dimensions(Letter_Mesh, FontWidth, FontLength, FontHeight);

		if (std::isdigit(c)) {
			lastWasDigit = true;
		}
		else if (lastWasDigit) {
			offsetY -= (FontLength * XYscale) + Yspacing;
			offsetX = -6.35; // 0.15
			lastWasDigit = false;
		}

		scaleMesh(Letter_Mesh, XYscale, XYtopscale, Zscale, zThreshold);
		translate_mesh(Letter_Mesh, Kernel::Vector_3(offsetX, offsetY, offsetZ + zDepth));
		offsetX += (FontWidth * XYscale) + Xspacing;
		CGAL::copy_face_graph(Letter_Mesh, Tag_Mesh);
	}

	Result_Mesh.clear();
	if (!PMP::corefine_and_compute_difference(Fixture_Mesh, Tag_Mesh, Result_Mesh)) {
		std::cerr << Red << "      Subtraction operation failed." << ColorEnd << std::endl;
	}
}

struct ModelType {
	std::string FullName;
	std::string label;
	int count;
};

bool processModel(const std::string outputPath, int ID, const ModelType modelType, int index) {
	std::string id = std::to_string(ID) + modelType.label + (index < 10 ? "0" : "") + std::to_string(index);
	std::string Filename = id + "_F.stl";
	std::string output = outputPath + "/" + Filename;

	std::cout << "      Creating: " << Filename << " for " << modelType.FullName << std::endl;

	Mesh Fixture_Mesh, Result_Mesh;

	if (!read_STL_data("fixture", Fixture_Mesh)) return false;

	create_fixture(id, Fixture_Mesh, Result_Mesh);

	if (!write_STL(output, Result_Mesh)) return false;
	return true;
}

void displayUserName() {
	char* username = nullptr;
	char* userdomain = nullptr;
	size_t sizeUsername = 0;
	size_t sizeUserdomain = 0;

	_dupenv_s(&username, &sizeUsername, "USERNAME");
	_dupenv_s(&userdomain, &sizeUserdomain, "USERDOMAIN");

	std::cout << "\n      USERNAME: "
		<< (userdomain ? userdomain : "Unknown") << "\\"
		<< (username ? username : "Unknown") << std::endl;

	// Free the allocated memory
	free(username);
	free(userdomain);
}

void promptForNumbers(const std::string& prompt, int& outValue) {
	std::string line;
	bool valid = false;

	while (!valid) {
		std::cout << Yellow << prompt << ColorEnd;
		std::getline(std::cin, line);
		if (line.empty()) {
			outValue = 0;
			valid = true;
		}
		else {
			std::stringstream ss(line);
			if (ss >> outValue && ss.eof()) {
				valid = true; 
			}
			else {
				std::cout << Red << "               Invalid input:  " << ColorEnd << line << std::endl;
				ss.clear(); 
			}
		}
	}
}


int main(int argc, char* argv[]) {

	//resizeConsole(77, 30);
	std::cout << Cyan << "=====================================================================" << ColorEnd << std::endl;
	std::cout << Cyan << "========================" << ColorEnd 
		<< Yellow << "'Created by Banna'" << ColorEnd
		<< Cyan << "===========================" << ColorEnd << std::endl;
	std::cout << Cyan << "===================" << ColorEnd 
		<< Yellow << "'AB FIXTURE CREATOR TOOL V3'" << ColorEnd
		<< Cyan<< "======================" << ColorEnd << std::endl;
	std::cout << Cyan << "=====================================================================\n" << ColorEnd << std::endl;


	std::string outputPath = fs::current_path().string() + "/output";
	if (!fs::exists(outputPath)) {
		fs::create_directory(outputPath);
	}
	else {
		for (const auto& entry : fs::directory_iterator(outputPath))
			fs::remove_all(entry.path());
	}

	int caseID;
	promptForNumbers("      What is the Case ID? ", caseID);

	std::vector<ModelType> models = {
		{"UPPER", "UN", 0}, {"UPPER PASSIVE", "UP", 0}, {"UPPER RETAINER", "UR", 0}, {"UPPER TEMPLATE", "UT", 0},
		{"LOWER", "LN", 0}, {"LOWER PASSIVE", "LP", 0}, {"LOWER RETAINER", "LR", 0}, {"LOWER TEMPLATE", "LT", 0}
	};

	for (auto& model : models) {
		promptForNumbers("      How many " + model.FullName + "? ", model.count);
	}

	auto start = std::chrono::high_resolution_clock::now();

	std::cout << Yellow << "\n============================'Creating Fixtures'==============================\n" << ColorEnd << std::endl;

	int processedCount = 0, totalCount = 0;
	for (const auto& model : models) {
		totalCount += model.count;
		for (int i = 1; i <= model.count; ++i) {
			if (processModel(outputPath, caseID, model, i)) {
				processedCount++;
			}
			else {
				std::cerr << Red << "      Failed to process " << ColorEnd
					<< model.FullName << " index " << i << std::endl;
			}
		}
	}
	std::cout << Yellow << "\n================================='Finished'==================================" << ColorEnd << std::endl;
	std::cout << Yellow <<   "=================================='REPORT'===================================\n" << ColorEnd << std::endl;



	std::cout << "      " << Green << processedCount << "/" << totalCount << ColorEnd
		<< "  OCR fixture created " << Green << "Successfully\n" << ColorEnd << std::endl;

	std::cout << "      " << Green << processedCount << ColorEnd 
		<< "  STL in 'output' - Fixture with OCR Tag\n" << std::endl;

	displayUserName();

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;

	std::cout << Yellow << "\n      Elapsed time: " << elapsed.count() << " seconds" << ColorEnd << std::endl;



	std::cout << "\nPress " << Green << "ENTER" << ColorEnd << " key to exit . . . " << std::endl;
	std::cin.get();
	return EXIT_SUCCESS;
}