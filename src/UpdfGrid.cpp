#include "UpdfGrid.h"
#include "Branchings.h"
#include <cmath>
#include "TH3.h"
#include "TFile.h"

#define UG_VERBOSE 1 
using namespace std;


// constructor
UpdfGrid::Updfgrid(double xmin, double xmax, double ktmin, double ktmax, double qmin, double qmax)
/* 
    parameters to be passed to make_grids  
*/
{
    make_grids(xmin, xmax, ktmin, ktmax, qmin, qmax);
}

// to build the grids
UpdfGrid::make_grids(double xmin, double xmax, double ktmin, double ktmax, double qmin, double qmax)
/* 
    (x/kt/q)min : lower bound for the (x/kt/scale) axis of the grid 
    (x/kt/q)max : upper bound for the (x/kt/scale) axis of the grid
*/
{
    double xbins[51] = {0};
    double ktbins[51] = {0};
    double qbins[51] = {0};
    // calculate the bin edges for the x axis of the grid
    inline void x_gridder(double (&bins)[51], double bmin, double bmax)
    {
      // define the special coarse spacing
      double xcoarse[6];
      xcoarse = {log(bmin), log(0.01), log(0.1), log(0.4), log(0.9), log(bmax)};
      // calculate the size of steps in the subdivisions
      double stepfine[5];
      for (int i=0; i<5; i++) {
        stepfine[i] = (xcoarse[i+1] - xcoarse[i])/10.;
      }
      // calculate the bin edges
      bins[0] = xcoarse[0];
      int incr = 0;
      for (int j=0; j<5; j++) {
        for (int k=0; k<10; k++) {
          bins[incr+1] = bins[incr] + stepfine[j];
          incr += 1;
        }
      }
    }
    // calculate the bin edges fot the kt or q axis of the grid
    inline void gridder(double (&bins)[51], double bmin, double bmax)
    {
      double step;
      double logbmin;
      logbmin = log(bmin);
      step = (log(bmax) - logbmin)/50.;
      for (int i=0; i<51; i++) {
        bins[i] = logbmin + step * i;
      }
    }

    // calculate the bin edges
    x_gridder(xbins, xmin, xmax);
    gridder(ktbins, ktmin, ktmax);
    gridder(qbins, qmin, qmax);

    // construct the 3d-grids
    *gluongrid = new("gluon_grid","gluons;x;kt;q",50,xbins,50,ktbins,50,qbins);
    *downgrid = new("down_grid","downs;x;kt;q",50,xbins,50,ktbins,50,qbins);
    *adowngrid = new("adown_grid","anti-downs;x;kt;q",50,xbins,50,ktbins,50,qbins);
    *upgrid = new("up_grid","ups;x;kt;q",50,xbins,50,ktbins,50,qbins);
    *aupgrid = new("aup_grid","anti-ups;x;kt;q",50,xbins,50,ktbins,50,qbins);
    *strangegrid = new("strange_grid","stranges;x;kt;q",50,xbins,50,ktbins,50,qbins);
    *astrangegrid = new("astrange_grid","anti-stranges;x;kt;q",50,xbins,50,ktbins,50,qbins);
    *charmgrid = new("charm_grid","charms;x;kt;q",50,xbins,50,ktbins,50,qbins);
    *acharmgrid = new("acharm_grid","anti-charms;x;kt;q",50,xbins,50,ktbins,50,qbins);
    *bottomgrid = new("bottom_grid","bottoms;x;kt;q",50,xbins,50,ktbins,50,qbins);
    *abottomgrid = new("abottom_grid","anti-bottoms;x;kt;q",50,xbins,50,ktbins,50,qbins);
    *topgrid = new("top_grid","anti-tops;x;kt;q",50,xbins,50,ktbins,50,qbins);
    *atopgrid = new("atop_grid","anti-tops;x;kt;q",50,xbins,50,ktbins,50,qbins);
    *photongrid = new("photon_grid","photons;x;kt;q",50,xbins,50,ktbins,50,qbins);
    *Zgrid = new("Z_grid","Zs;x;kt;q",50,xbins,50,ktbins,50,qbins);
    *Wminusgrid = new("Wminus_grid","Wminuss;x;kt;q",50,xbins,50,ktbins,50,qbins);
    *Wplusgrid = new("Wplus_grid","Wpluss;x;kt;q",50,xbins,50,ktbins,50,qbins);
    *Higgsgrid = new("Higgs_grid","Higgss;x;kt;q",50,xbins,50,ktbins,50,qbins);
}

// to fill the grids
UpdfGrid::add_event(parton event)
/* 
    event : parton whose data is to be stored in one of the grids 
*/
{
    double ktx, kty, x, kt2, weight, scale;
    int flavour;

    // get the data from the event
    ktx = event.ktx;
    kty = event.kty;
    x = event.x;
    kt2 = ktx*ktx + kty+kty;
    weight = event.weight;
    q = event.scale;
    flavour = event.flavour;

    // fill the grid of the corresponding flavour
    switch (flavour)
    {
    case 0:
      gluongrid->Fill(log(x),log(kt2),log(q),weight/x); break;
    case 1:
      downgrid->Fill(log(x),log(kt2),log(q),weight/x); break;
    case -1:
      adowngrid->Fill(log(x),log(kt2),log(q),weight/x); break;
    case 2:
      upgrid->Fill(log(x),log(kt2),log(q),weight/x); break;
    case -2:
      aupgrid->Fill(log(x),log(kt2),log(q),weight/x); break;
    case 3:
      strangegrid->Fill(log(x),log(kt2),log(q),weight/x); break;
    case -3:
      astrangegrid->Fill(log(x),log(kt2),log(q),weight/x); break;
    case 4:
      charmgrid->Fill(log(x),log(kt2),log(q),weight/x); break;
    case -4:
      acharmgrid->Fill(log(x),log(kt2),log(q),weight/x); break;
    case 5:
      bottomgrid->Fill(log(x),log(kt2),log(q),weight/x); break;
    case -5:
      abottomgrid->Fill(log(x),log(kt2),log(q),weight/x); break;
    case 6:
      topgrid->Fill(log(x),log(kt2),log(q),weight/x); break;
    case -6:
      atopgrid->Fill(log(x),log(kt2),log(q),weight/x); break;
    case 7:
      photongrid->Fill(log(x),log(kt2),log(q),weight/x); break;
    case 8:
      Zgrid->Fill(log(x),log(kt2),log(q),weight/x); break;
    case 9:
      Wminusgrid->Fill(log(x),log(kt2),log(q),weight/x); break;
    case 10:
      Wplusgrid->Fill(log(x),log(kt2),log(q),weight/x); break;
    case 11:
      Higgsgrid->Fill(log(x),log(kt2),log(q),weight/x); break;
    }  
}