#ifndef VTK_EXPORTER_H
#define VTK_EXPORTER_H

#include <string>
#include <eigen3/Eigen/Core>

// Exportador VTK para campos escalares definidos sobre una malla
// estructurada regular (p. ej. la temperatura de una placa 2D).
//
// A diferencia de una simulación de partículas (DATASET UNSTRUCTURED_GRID,
// un punto por partícula), aquí el estado es un campo definido en los
// nodos de una malla regular Nx * Ny, por lo que el formato natural de
// VTK es STRUCTURED_POINTS: basta con indicar dimensiones, origen y
// espaciado, y listar los valores del campo en orden (x variando más
// rápido, luego y)

class VTKExporter {
private:
    std::string output_dir;
    int frame_counter = 0;

    int nx;      // Número de nodos en X
    int ny;      // Número de nodos en Y
    double dx;   // Espaciado entre nodos en X
    double dy;   // Espaciado entre nodos en Y

public:
    // Constructor: directorio de salida y geometría de la malla
    explicit VTKExporter(const std::string& output_directory,
                          int nx, int ny, double dx, double dy);

    // Exporta el campo escalar (p. ej. temperatura) del instante actual
    // como un archivo .vtk legible por ParaView.
    void write_frame(const Eigen::VectorXd& field, double time,
                      const std::string& field_name = "temperatura");

    // Reinicia el contador de frames (por si se corre otra simulación)
    void reset_counter();
};

#endif
