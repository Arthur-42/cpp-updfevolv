#ifndef UPDFGRID_H
#define UPDFGRID_H
#include "TH3.h"
#include "TFile.h"
#include "Branchings.h"

class UpdfGrid
{
    private:
        /* grids (as 3d-histograms) */ //possiblement à compléter
        TH3D *gluongrid;
        TH3D *downgrid;
        TH3D *adowngrid;
        TH3D *upgrid;
        TH3D *aupgrid;
        TH3D *strangegrid;
        TH3D *astrangegrid;
        TH3D *charmgrid;
        TH3D *acharmgrid;
        TH3D *bottomgrid;
        TH3D *abottomgrid;
        TH3D *topgrid;
        TH3D *atopgrid;
        TH3D *photongrid;
        TH3D *Zgrid;
        TH3D *Wminusgrid;
        TH3D *Wplusgrid;
        TH3D *Higgsgrid;

    public:
        /* methods */
        // constructor
        Updfgrid(double xmin, double xmax, double ktmin, double ktmax, double qmin, double qmax);
        // to build the grids
        void make_grids(double xmin, double xmax, double ktmin, double ktmax, double qmin, double qmax);
        // to fill the grids
        void add_event(parton event);
        // to save the grids to a file
        void save_grids();
}

#endif