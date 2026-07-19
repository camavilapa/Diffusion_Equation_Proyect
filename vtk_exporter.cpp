#include "vtk_exporter.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

VTKExporter::VTKExporter(const std::string& output_directory,
                          int nx_, int ny_, double dx_, double dy_)
    : output_dir(output_directory), nx(nx_), ny(ny_), dx(dx_), dy(dy_) {
    // Crea el directorio de salida si no existe
    if (!fs::exists(output_dir)) {
        try {
            fs::create_directories(output_dir);
            std::cout << "Directorio de salida creado: " << output_dir << "\n";
        } catch (const std::exception& e) {
            std::cerr << "Advertencia: no se pudo crear el directorio de salida: "
                      << e.what() << "\n";
        }
    }
}

void VTKExporter::write_frame(const Eigen::VectorXd& field, double time,
                               const std::string& field_name) {
    // Nombre del archivo: frame_000000.vtk, frame_000001.vtk, ...
    std::ostringstream filename_stream;
    filename_stream << output_dir << "/frame_"
                     << std::setfill('0') << std::setw(6) << frame_counter
                     << ".vtk";
    std::string filename = filename_stream.str();

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: no se pudo abrir el archivo para escritura: "
                   << filename << "\n";
        return;
    }

    // --- Encabezado formato VTK Legacy ---
    file << "# vtk DataFile Version 2.0\n";
    file << "Difusion de calor 2D - Tiempo: " << std::fixed << std::setprecision(6)
         << time << " s\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";

    // Malla regular: nx * ny nodos, 1 sola capa en Z (problema 2D)
    file << "DIMENSIONS " << nx << " " << ny << " 1\n";
    file << "ORIGIN 0.0 0.0 0.0\n";
    file << "SPACING " << dx << " " << dy << " 1.0\n";

    // --- Campo escalar en cada nodo (p. ej. temperatura) ---
    // El orden esperado por VTK es x variando más rápido y luego y,
    // que coincide con k = i*nx + j usado en pruebas.cpp.
    file << "POINT_DATA " << (nx * ny) << "\n";
    file << "SCALARS " << field_name << " float 1\n";
    file << "LOOKUP_TABLE default\n";

    for (int k = 0; k < nx * ny; ++k) {
        file << std::scientific << std::setprecision(8) << field[k] << "\n";
    }

    file.close();
    std::cout << "Frame VTK escrito: " << filename << " (t = " << time << " s)\n";
    frame_counter++;
}

void VTKExporter::reset_counter() {
    frame_counter = 0;
}
