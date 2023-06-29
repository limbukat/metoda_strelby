#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

using std::cin;
using std::cout;
using std::vector;

bool strelba(double z_0, double y_m, double y0_d, double y0_h, double epsilon, double deltaX);

const double xMin = 0;
const double xMax = 1;
int pocetBodu;

int main()
{
    double alpha, beta;
    double z0_dolni, z0_horni;
    double deltaX;
    double epsilon = 1e-6;

    cout << "NME: Uloha 4\n";
    cout << "y(0) = alpha\n";
    cout << "y(1) = beta\n";
    cout << "epsilon = 1e-6\n\n";

    cout << "Zadej alpha: ";
    cin >> alpha;
    cout << "Zadej beta: ";
    cin >> beta;
    cout << "Zadej dolni odhad y'(0): ";
    cin >> z0_dolni;
    cout << "Zadej horni odhad y'(0): ";
    cin >> z0_horni;
    cout << "Zadej delku kroku: ";
    cin >> deltaX;

    pocetBodu = ((xMax-xMin)/deltaX) +1;
    strelba(alpha, beta, z0_dolni, z0_horni, epsilon, deltaX);
}

double zDerivace(double x, double y) // zDerivace = yDruhaDerivace (protoze z := y'), y''=f(x,y,y')
{
    return ((-21*x*x*x*x - 36*x*x*x*x*x*x*x + 36*x*x*x*x*x*x*x*x*x*x)*y + 6*x*pow(2.718281828,x*x*x*(1-x*x*x)));
}

vector<double> euler( double y_0, double z_0, double deltaX)   // Eulerova metoda - øešení dif rovnice 
{
    vector<double> x(pocetBodu);
    vector<double> y(pocetBodu);
    vector<double> z(pocetBodu);

    x[0] = xMin;
    y[0] = y_0;
    z[0] = z_0;

    for (int i=0; i < pocetBodu-1; i++)
    {
        x[i+1] = x[0] + deltaX * (i+1);
        y[i+1] = y[i] + deltaX * z[i];
        z[i+1] = z[i] + deltaX * zDerivace(x[i],y[i]);
    }
    return y;
}

bool jsouPocatecniOdhadyDobreDefinovane(vector<double> y_dolni, vector<double> y_horni, double beta)   // využiju pøi metodì pùlení intervalù
{
    if( y_dolni[pocetBodu-1] > beta )
    {
        cout << "Pro dolni odhad je y(1) vetsi nez beta\n";
        return false;
    }
    if ( y_horni[pocetBodu-1] < beta )
    {
        cout << "Pro horni odhad je y(1) mensi nez beta\n";
        return false;
    }
    return true;
}

bool jeSpoctenaFunkceSpravnymResenim(vector<double> y, double beta, double epsilon)   // kdy ukonèím metodu pùlení intervalu, strefila jsem se? 
{
    return (std::abs(y[pocetBodu-1] -beta) < epsilon);
}

bool strelba(double y_0, double beta, double z0_dolni, double z0_horni, double epsilon, double deltaX) //metoda støelby 
{
    vector<double> x(pocetBodu);
    vector<double> y_dolni;
    vector<double> y_horni;

    for (int i=0; i < pocetBodu; i++)    // vytvoøím vektor bodù x 
        x[i]=i*deltaX;

    std::ofstream fout("NME_uloha4.txt");
    int iterace = 0;
    while (true)
    {
        ++iterace;
        y_dolni = euler(y_0, z0_dolni, deltaX);    // y(b, odhad1)
        y_horni = euler(y_0, z0_horni, deltaX);    // y(b, odhad2)

        if (!jsouPocatecniOdhadyDobreDefinovane(y_dolni, y_horni, beta))   // musím zadat nový odhad
            return false;

        if (jeSpoctenaFunkceSpravnymResenim(y_dolni, beta, epsilon))  // trefila jsem se?  
        {
            for (int i=0; i<pocetBodu; i++)
                fout << x[i] << "\t" << y_dolni[i] << "\n";
            return true;
        }

        if (jeSpoctenaFunkceSpravnymResenim(y_horni, beta, epsilon))  // trefila jsem se?
        {
            for (int i=0; i<pocetBodu; i++)
                fout << x[i] << "\t" << y_horni[i] << "\n";
            return true;
        }

        double z0_novyOdhad = (z0_dolni + z0_horni)/2;  // alfa*
        vector<double> y_novyOdhad = euler(y_0, z0_novyOdhad , deltaX);   // y(b,alfa*)
        if (jeSpoctenaFunkceSpravnymResenim(y_novyOdhad, beta, epsilon))   // trefila jsem se? 
        {
            for (int i=0; i<pocetBodu; i++)
                fout << x[i] << "\t" << y_novyOdhad[i] << "\n";
            return true;
        }

        if (y_novyOdhad[pocetBodu-1] > beta)
            z0_horni = z0_novyOdhad;
        else
            z0_dolni = z0_novyOdhad;
    }
}
