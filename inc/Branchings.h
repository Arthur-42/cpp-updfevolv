#ifndef BRANCHINGS_H
#define BRANCHINGS_H

// holds the relevant data at each step of the evolution
struct parton
{
    double ktx;      // kt_x of the parton
    double kty;      // kt_y of the parton
    double x;          // current momentum fraction
    double weight;     // current full weight
    double scale;      // current scale
    int flavour;   /*indicates the type of the particle : */
};

// manages the evolution of one parton from starting to final scale
class Branchings
{
    private:
        /* data */
        // to hold the data of the current branching
        parton m_parton;
        // internal variables
        double m_xmin;
        double m_xmax;
        double m_starting_scale;

        /* methods */
        // to generate the starting parton
        void start_the_chain(double (*starting_pdf)(double), int starting_flavour, 
                double kt0_min, double kt0_sigma);    
        // to generate the scale of the next branching
        void generate_next_scale();    
        // ...

    public:
        /* methods */
        // construtor
        Branchings(double (*starting_pdf)(double), double starting_scale, int starting_flavour, 
                double xmin, double xmax, double kt0_min, double kt0_sigma);    
        // to manage the evolution process              
        void evolution(&grids);        
};

#endif