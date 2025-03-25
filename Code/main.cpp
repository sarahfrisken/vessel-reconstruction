#include "VesselReconstruction.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    VesselReconstruction w;
    w.show();
    return a.exec();
}
