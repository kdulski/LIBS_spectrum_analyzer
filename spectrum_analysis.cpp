#include <iostream>
#include <fstream>
#include <string>
#include "TH1F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "THStack.h"
#include "TMath.h"
#include "TF1.h"
#include "TFile.h"

using namespace std;

Double_t start_ch1 = 182.56;
Double_t stop_ch1 = 252.94;
Double_t start_ch2 = 252.95;
Double_t stop_ch2 = 314.30;
Double_t start_ch3 = 314.31;
Double_t stop_ch3 = 419.17;
Double_t start_ch4 = 419.18;
Double_t stop_ch4 = 498.42;
Double_t start_ch5 = 498.43;
Double_t stop_ch5 = 721.85;
Double_t start_ch6 = 721.86;
Double_t stop_ch6 = 909.27;

double tempMeanSD = 1;

int minLength_ch[6] = {4, 4, 4, 6, 4, 4}; // ile binów musi być spadek w każdym z kanałów
int halfLengthToCheckExtremum_ch[6] = {5, 5, 5, 5, 3, 3}; // ile punktow obok musi byc mniejsze od maxa

TH1F* Prepare_spectrum(vector<double> arguments, vector<double> values, double start, double stop, int ch);
double GramPolynomial(int i, int m, int k, int s);
double GenFactor(int a, int b);
double CalcWeight(int i, int t, int m, int n, int s);
vector<double> ComputeWeights(int sizeOfHalfWindow, int numberOfPoints, int orderOfPolynomialToInterpolate, int degreeOfDerivative);
TGraph* FindDecreasingSeries(TH1F* derivative_hist, int minLength, double halfLengthToCheckExtremum, int minimumHeightDifference);
double VoigtForPeak(double x, double mean, double sigma, double gamma);
double VoigtToDraw(Double_t *A, Double_t *P);
double StandardDeviation(TGraph* spectrum_ch);
TGraph* EstimateBackground(vector<double> arguments_ch, vector<double> values_ch);
TGraph* FindPeaks(vector<double> arguments_ch, vector<double> values_ch, TF1* bgFit);
Double_t peakFunc(Double_t *A, Double_t *P);

vector<double> peakWavelenght; // globalny wektor na dł. fali peaków
vector<double> peakWidth; // globalny wektor na szerokości peaków
vector<double> peakHeight;
// double InitialNoiseRangeNm = 190; // będzie liczył odchylenie standardowe z zakresu x od początku do x=InitialNoiseRange, tu jest dł. fali

int main (int argc, char* argv[]) {
  string data_file;
  if (argc < 2) {
    cout << "No arguments given. Turning on interactive input mode." << endl;
    cout << "Give me the name of the file to analyze" << endl;
    cin >> data_file;
  } else {
    data_file = argv[1];
  }

// deklaracja obiektu typu fstream
    fstream file;
// otwarcie pliku w trybie ios::in czyli do odczytu
    file.open(data_file, ios::in);

    if(!file.is_open()){
        cout << "Error while opening file: " << data_file << endl;
        return 1;
    } else {
        cout << "File: " << data_file << " opened" << endl;
    }

    string line;
    string tempArg, tempValue;
// tworzy 2 wektory domyślnie puste
    vector<double> arguments, values;

// wczytuje wiersz z "file" i zapisuje go jako "line"
    getline(file, line);
// założenie że mamy dwie liczby na wiersz gdzie pierwsza wartosc z wiersza jest zapisywana jako tempArg a druga jako tempValue
// jeśli w pliiku jest jakiś nagłówek trzeba ominąć te linie żeby ominąc występowanie błędu
    while (file >> tempArg >> tempValue) {
// dodaje wartosc z stod(tempArg) na koniec wektora arguments
        arguments.push_back(stod(tempArg));
        values.push_back(stod(tempValue));
    }

    vector<pair<double,double>> limits_ch;
    limits_ch.push_back(make_pair(start_ch1, stop_ch1));
    limits_ch.push_back(make_pair(start_ch2, stop_ch2));
    limits_ch.push_back(make_pair(start_ch3, stop_ch3));
    limits_ch.push_back(make_pair(start_ch4, stop_ch4));
    limits_ch.push_back(make_pair(start_ch5, stop_ch5));
    limits_ch.push_back(make_pair(start_ch6, stop_ch6));

    TGraph *spectrum[6];
    TCanvas* c1[6];

    double sigma[6];
    double minimumHeightDifference_ch[6];
// wektory na argumenty i wartości dla kanału w którym jesteśmy
    vector<double> arguments_ch, values_ch;

    double noOfSigmaForConfInterval = 1;

    TString slash = "/";
    TString dot = ".";
    std::size_t slashPlace = data_file.rfind(slash);
    std::size_t dotPlace = data_file.rfind(dot);
    std::string filenameCleaned = data_file;
    if (slashPlace > 0 && dotPlace > slashPlace) {
      filenameCleaned = data_file.substr(slashPlace+1, dotPlace);
    } else if (dotPlace > 0) {
      filenameCleaned = data_file.substr(0, dotPlace);
    }
    TString outputName = "LIBSFit_" + filenameCleaned + ".root";
    TFile* outputFile = new TFile(outputName, "RECREATE");

    for (Int_t i=0; i<6; i++) {
        arguments_ch.clear();
        values_ch.clear();

        for (Int_t j=0; j<arguments.size(); j++) { // dla kolejnych argumentów (wszystkich)
            if (arguments[j] >= limits_ch.at(i).first && arguments[j] <= limits_ch.at(i).second) {
// jeżeli argument należy do kanału w którym jesteśmy
                arguments_ch.push_back(arguments[j]);
                values_ch.push_back(values[j]);
            }
        }

// spectrum[i] = new TH1F(Form("Spectrum channel %d", i + 1), Form("Spectrum channel %d", i + 1), nbins, start_ch[i], stop_ch[i]); // (liczba binow - przedzialow, dolna krawedz pierwszego binu, gorna krawedz ostatniego binu)
        TGraph *spectrum = new TGraph(arguments_ch.size(), arguments_ch.data(), values_ch.data());
        spectrum->Write("rawSpectrum");

// (tytul w root, tytuł w oknie, szerokość w pikselach, wysokość w pikselach)
//        c1[i] = new TCanvas(Form("c1 ch%d", i + 1), Form("c1 ch%d", i + 1), 800, 600);

//        spectrum[i]->SetLineColor(kBlue+3);
//        spectrum[i]->SetLineWidth(2);
// SetOptStat okresla jakie informacje  maja sie wyswietlic na wykresie
//        gStyle->SetOptStat(111111);

        // spectrum[i]->Draw("HIST"); // rysuj słupki
//        spectrum[i]->Draw("AL"); // A - axes, L -lines

//        sigma[i] = StandardDeviation(spectrum[i]);
        TGraph *bgGraph = EstimateBackground(arguments_ch, values_ch);
        TF1* bgFit = new TF1("bgFit", "pol4", (Double_t)limits_ch.at(i).first, (Double_t)limits_ch.at(i).second);
        bgGraph->Fit("bgFit");
        bgGraph->Write("bgGraph");

        TGraph *peaksCandidatesGraph = FindPeaks(arguments_ch, values_ch, bgFit);
        peaksCandidatesGraph->Write("peakGraph");

        int noOfPeaks = peaksCandidatesGraph->GetN();
        if (noOfPeaks > 0) {
          TF1* peakFit = new TF1("peakFit", peakFunc, (Double_t)limits_ch.at(i).first, (Double_t)limits_ch.at(i).second, 8+2*noOfPeaks);
          peakFit->SetParameter(0, bgFit->GetParameter(0));
          peakFit->SetParameter(1, bgFit->GetParameter(1));
          peakFit->SetParameter(2, bgFit->GetParameter(2));
          peakFit->SetParameter(3, bgFit->GetParameter(3));
          peakFit->SetParameter(4, bgFit->GetParameter(4));
          peakFit->FixParameter(5, noOfPeaks);
          peakFit->SetParameter(6, 0.1); // gamma
          peakFit->SetParameter(7, 0.1); // sigma
          double LargeIntens = 1E10;
          for (unsigned j=0; j<noOfPeaks; j++) {
            peakFit->SetParameter(8+j*2, peaksCandidatesGraph->GetPointY(j));
            peakFit->SetParName(8+j*2, "Intensity");
            peakFit->SetParLimits(8+j*2, 0, LargeIntens);

            peakFit->SetParameter(9+j*2, peaksCandidatesGraph->GetPointX(j));
            peakFit->SetParName(9+j*2, "Mean");
            peakFit->SetParLimits(9+j*2, (Double_t)limits_ch.at(i).first, (Double_t)limits_ch.at(i).second);
          }
          cout << "Fitting " << noOfPeaks << " peak(s)" << endl;
          spectrum->Fit(peakFit, "RM");
          peakFit->SetNpx(1000000);
          TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
          spectrum->Draw("");
          peakFit->Draw("same");
          c1->Write("FittedSpectrum");
          delete c1;
        }
//        cout << "sigma ch " << i+1 << " = " << sigma[i] << endl;

//        minimumHeightDifference_ch[i] = noOfSigmaForConfInterval * sigma[i]; // różnica w wysokości max i min pochodnej
    }

    outputFile->Close();
    return 0;

    TH1F* derivative[6];
    TGraph* peaks[6];
    THStack* stack = new THStack("stack", "Smoothed derivatives; wavelength [nm]; counts [a.u.]");

    for (Int_t i=0; i<6; i++) {
        derivative[i] = Prepare_spectrum(arguments, values, limits_ch.at(i).first, limits_ch.at(i).second, i+1); // wywolanie funkcji
        derivative[i]->SetTitle(Form("Smoothed devivative - channel no. %d", i + 1)); // tytuł w oknie wykresu
        derivative[i]->SetLineColor(kRed+3);
        derivative[i]->SetLineWidth(4);
        gStyle->SetOptStat(111111); // SetOptStat okresla jakie informacje  maja sie wyswietlic na wykresie

        stack->Add(derivative[i]);

        peaks[i] = FindDecreasingSeries(derivative[i], minLength_ch[i], halfLengthToCheckExtremum_ch[i], minimumHeightDifference_ch[i]);

        cout << "Liczba znalezionych peaków w kanale " << i+1 << " to " << peaks[i]->GetN() << endl; // GetN() - liczba punktow w graphie

        peaks[i]->SetMarkerStyle(8); // to kropka
        peaks[i]->SetMarkerColor(kRed);
        peaks[i]->SetMarkerSize(0.5);

// rysowanie peaków na spektrum
        c1[i]->cd();
        peaks[i]->Draw("P SAME");
    }

    TCanvas *c2 = new TCanvas("c2", "Smoothed Derivative", 800, 600);

    stack->Draw("NOSTACK HIST");   // NOSTACK HIST - wartości z histogramów się nie sumują

    stack->SetMinimum(-40000); // ustawienie zakresu osi y
    stack->SetMaximum(40000);

    TCanvas *c3 = new TCanvas("c3", "Peak Widths", 800, 600); 

    cout << "peakWavelenght = " << peakWavelenght.size() << endl;
    // cout << "peakWidth = " << peakWidth.size() << endl;

    TGraph *peakwidths = new TGraph(
    peakWavelenght.size(),
    peakWavelenght.data(),
    peakWidth.data()
    );

    peakwidths->SetTitle("Peakwidths;peak wavelength [nm];peak width [nm]");
    peakwidths->SetMarkerStyle(20);
    peakwidths->SetMarkerSize(1.0);
    peakwidths->SetMarkerColor(kBlue+3);

    peakwidths->Draw("AP"); // osie i punkty


    
    for (int j = 0; j < 6; j++) {
        c1[j]->cd();
        for (Int_t i = 0; i < peakWavelenght.size(); i++) { // pętla po peakach

            if (peakWavelenght[i] >= limits_ch.at(j).first &&
            peakWavelenght[i] <= limits_ch.at(j).second) {

                TF1 *voigtFit = new TF1(Form("voigtFit_ch%d_peak%d", j + 1, i), VoigtToDraw, peakWavelenght[i] - peakWidth[i], peakWavelenght[i] + peakWidth[i], 6); // (nazwa, wzór, odkąd, dokąd, ile parametrów)

                voigtFit->SetParameters(
                peakHeight[i], // amplituda
                peakWavelenght[i], // mean, czyli x w środku
                peakWidth[i]/2,   // sigma, 
                peakWidth[i]/4,  // lg
                1,
                1
                );
// dodać sprawdzenie które spektrum - element w tablicy/wektorze; na podstawie start i stop ch
// for int j=0 do 6
//    if (peakWavelenght[i] > startCh[j] && peakWavelenght[i] < stop_ch[j])
                spectrum[j]->Fit(voigtFit,"RM+"); // R - fituj w ustalonym regionie, + - nie usuwaj poprzednich fitów
                voigtFit->Draw("SAME");
            }
        }
    }
    
    return 0;
}


TH1F* Prepare_spectrum(vector<double> arguments, vector<double> values, double start, double stop, int ch) // histogram z wygladzona pochodna spektrum
{
    int sizeOfHalfWindow = 5; // nie isc wyzej niz 8
    int orderOfPolynomialToInterpolate = 2; // mozna zmienic na 3

    vector<double> weightsForDerivative = ComputeWeights(sizeOfHalfWindow, 0, orderOfPolynomialToInterpolate, 1);

    int iterator = 1;
    while (arguments.at(iterator) < start) {
        iterator++;
    }
    double binSize = arguments.at(iterator+1) - arguments.at(iterator);
    TH1F* derivative_hist = new TH1F(Form("derivative_hist_%d", ch),
    Form("derivative_hist_%d", ch), (stop-start)/binSize, start, stop);

    for(unsigned i = iterator + sizeOfHalfWindow; i < arguments.size()-sizeOfHalfWindow && arguments.at(i+sizeOfHalfWindow) < stop; i++){ // idziemy po kolejnych punktach na osi x od sizeOfHalfWindow
        double tempValueOfCounts = 0; // tymczasowa zmienna dla kolejnych iloczynów counts * waga
        for(int j = -sizeOfHalfWindow; j <= sizeOfHalfWindow; j++){
            tempValueOfCounts += weightsForDerivative.at(j+sizeOfHalfWindow)*values.at(i+j);
        }
        int bin = derivative_hist->FindBin(arguments[i]); // FindBin(arguments[i]) znajduje bin, w którym jest arguments[i]
        derivative_hist->SetBinContent(bin, tempValueOfCounts); // wypelniamy bin o nazwie "bin" wartoscią po przecinku
    }

    return derivative_hist;
}



double GramPolynomial(int i, int m, int k, int s)
{
    if (k>0) {
        return (4.*k - 2.)/(k*(2.*m - k + 1.))*(i*GramPolynomial(i, m, k-1, s) + s*GramPolynomial(i, m, k-1, s-1)) - ((k-1.)*(2.*m + k))/(k*(2.*m - k + 1.))*GramPolynomial(i, m, k-2, s);
    } else {
        if (k==0 && s==0)
            return 1.;
        else
            return 0.;
    }
}

double GenFactor(int a, int b)
{
    double gf = 1.;
    for (int j=(a-b)+1; j<=a; j++) {
        gf*=j;
    }
    return gf;
}

double CalcWeight(int i, int t, int m, int n, int s)
{
    double w = 0;
    for (int k=0; k<=n; ++k) {
        w = w + (2*k + 1)*(GenFactor(2*m, k)/GenFactor(2*m + k + 1, k + 1))*GramPolynomial(i, m, k, 0)*GramPolynomial(t, m, k, s);
    }
    return w;
}

// m -> size of half window - how many points from the central are taken
// t -> number of data point (give 0)
// n -> order of interpolation polynomial -> 1 - line, 2 - quadratic ...
// s -> order of function, 0 - smooth, 1 - derivative, 2 - second derivative, ...
std::vector<double> ComputeWeights(int m, int t, int n, int s)
{
    std::vector<double> weights(2*m + 1);
    for (int i=0; i<2*m+1; ++i) {
        weights.at(i) = CalcWeight(i-m, t, m, n, s);
    }

  //  std::cout << "Calculated weights ";
 //   for (unsigned i=0; i<weights.size(); i++)
   //     std::cout << weights.at(i) << " ";
   // std::cout << std::endl;

    return weights;
}

TGraph* FindDecreasingSeries(TH1F* derivative_hist, int minLength, double halfLengthToCheckExtremum, int minimumHeightDifference)
{
    // importance of maximum will be defined by the length and the amplitude of the constantly decreasing series
    int startBin = 0, stopBin = 0; // startBin to ostatnie znalezione lokalne maksimum pochodnej, stopBin to ostatnie znalenione min
    // int minLength = 3; // zmniejszyłam z 10 - ile binow musi byc spadek
    // int halfLengthToCheckExtremum = 5; // zmniejsztlam z 11 - ile punktow obok musi byc mniejsze od maxa
    bool testMax, testMin;
    double extremaHeightDifference; // to zmienna do przechowywania różnicy pomiędzy wysokością ostatnio znaleznionego max i min
    // int MinimumHeightDifference = 1800; // min. różnica w wysokości max i min pochodnej

    std::vector<std::pair<int,int>> localExtrema; // wektor do gromadzenia par: ostatnie max pochodnej + ostatnie min pochodnej

    for (unsigned i=halfLengthToCheckExtremum+1; i<derivative_hist->GetNbinsX() - halfLengthToCheckExtremum; i++) { 
        if (startBin != 0) {
            if (derivative_hist->GetBinContent(startBin) < derivative_hist->GetBinContent(i)) // jezeli wartosc pochodnej w i jest większa niż w startBin -> jest wzrost-> ustawiamy startBin na 0
                startBin = 0;
        }

        if (startBin == 0) {
            testMax = true;
            int j=1; // j jest po to żeby iść po sąsiednich binach
            while (testMax && j<=halfLengthToCheckExtremum) { // jeżeli testMax jest true (czyli startBin = 0) i j jest w granicach w jakich sprawdzamy spadek
                if ((derivative_hist->GetBinContent(i) < derivative_hist->GetBinContent(i+j) || derivative_hist->GetBinContent(i) < derivative_hist->GetBinContent(i-j)) && testMax) { // czy wartosc w którymkolwiek z sąsiednich binów jest większa niż w i
                    testMax = false; // jeżeli tak jest to nie jesteśmy w maximum więc testMax - false
                }
                j++; // sprawdzamy kolejnych sąsiadów
            }
            // jeżeli punkt i przejdzie ten test to uznajemy go za max pochodnej

            if (testMax) {
                startBin = i; // jezeli punkty dookola nie sa mniejsze to ustawiamy max ("startBin") na i <- znaleźliśmy peak pochodnej

            }
        }
        else { // jeżeli startBin =/= 0, czyli nie ma spadku od ostatniego maxa pochodnej
            testMin = true; 
            int j=1;
         /*   while (testMin && j<=halfLengthToCheckExtremum) {
                if ((derivative_hist->GetBinContent(i) > derivative_hist->GetBinContent(i+j) || derivative_hist->GetBinContent(i) > derivative_hist->GetBinContent(i-j)) && testMin)
                    testMin = false;
                j++;
            }
            */
            if (i < derivative_hist->GetNbinsX() - halfLengthToCheckExtremum){ // czy nie wyszliśmy poza wykres
                while(testMin && j<=halfLengthToCheckExtremum){ 
                    if(derivative_hist->GetBinContent(i) > derivative_hist->GetBinContent(i-j) || derivative_hist->GetBinContent(i) > derivative_hist->GetBinContent(i+j)) { // czy sasiednie sa wieksze
                        testMin = false;
                    }
                    else{
                        testMin = true;
                        j++; // sprawdzamy kolejnych sasiadów
                    }
                }
            }
            
            if (testMin) {
                stopBin = i; // jezeli punkty dookola nie sa mniejsze to ustawiamy min ("stopBin") na i <- znaleźliśmy min pochodnej

                extremaHeightDifference =  derivative_hist->GetBinContent(startBin) - derivative_hist->GetBinContent(stopBin); // liczymy różnicę w wysokości max i min pochodnej

                if (i - startBin > minLength && extremaHeightDifference >= minimumHeightDifference) { // czy odległość od min pochodnej do max pochodnej jest wystarczajaco duża
                    localExtrema.push_back(std::make_pair(startBin, stopBin)); // wpisujemy wartości startBin i stopBin na koniec wektora localExtrema
                    cout << "Znalezniono max pochodnej dla x = " << startBin
                    << " i min pochodnej dla x =" << stopBin << endl;
                }
                startBin = 0;
            }
            
            //cout << derivative_hist->GetBinCenter(startBin) << " " << derivative_hist->GetBinCenter(i) << endl;
           /* if (!testMin) {
                if (i - startBin > minLength) { // czy odległość od min pochodnej do max pochodnej jest wystarczajaco duża
                    localMaxima.push_back(std::make_pair(startBin, stopBin));
                }
                startBin = 0;
            }
            */    
        }
    }
    const Int_t n = localExtrema.size();
    vector<Double_t> x(n), y(n); // zmieniłam tablice na wektory bo było ostrzeżenie
    for (unsigned i=0; i<localExtrema.size(); i++) {
        //int centralBin  = (localExtrema.at(i).second + localExtrema.at(i).first)/2; // środek prostej między max i min pochodnej
        // y = ax + b
        double a = (derivative_hist->GetBinContent(localExtrema.at(i).first) - derivative_hist->GetBinContent(localExtrema.at(i).second)) / // first - pierwsza wartość z pary, second - druga wartość z pary
        (derivative_hist->GetBinCenter(localExtrema.at(i).first) - derivative_hist->GetBinCenter(localExtrema.at(i).second));
        double b = derivative_hist->GetBinContent(localExtrema.at(i).first) - a*derivative_hist->GetBinCenter(localExtrema.at(i).first);
        x[i] = -b/a;
        y[i] = - 1000;
        // y[i] = derivative_hist->GetBinContent(localMaxima.at(i).first) - derivative_hist->GetBinContent(localMaxima.at(i).second);
        cout << "Znaleziono peak dla x = " << x[i] << endl;
        peakWavelenght.push_back(x[i]);
        peakHeight.push_back(derivative_hist->GetBinContent(localExtrema.at(i).first)); // dodaje startBin do wektora wysokości peaków

        // Liczenie szerokości peaku i dodawanie tej wartości do wektora
        double width = derivative_hist->GetBinCenter(localExtrema.at(i).second) - derivative_hist->GetBinCenter(localExtrema.at(i).first);
        peakWidth.push_back(width);
    }

    return new TGraph(n,x.data(), y.data()); // zmieniłam x i y na x.data() i y.data() czyli wskaźniki do danych z wektorów
}


double VoigtForPeak(double x, double mean, double sigma, double gamma)
{
  return sqrt(2*M_PI)*sigma*TMath::Voigt(x-mean, sigma, gamma, 4);
}


double VoigtToDraw(Double_t *A, Double_t *P)
{
  double amplitude = P[0];
  double energy = P[1];
  double sigma = P[2];
  double gamma = P[3];

  double aParam = P[4];
  double bParam = P[5];

  double background = aParam*A[0] + bParam;
  return background + amplitude*VoigtForPeak(A[0], energy, sigma, gamma);
//  return P[0]*sqrt(2*M_PI)*sigma*TMath::Voigt(A[0] - energy, sigma, P[3], 4);
}


double StandardDeviation(TGraph* spectrum_ch)
{   
    // int InitialNoiseRangeNm = 20
    // int InitialNoiseRange = hist->FindBin(InitialNoiseRangeNm); // liczy jaki to bin z dł. fali podanej na początku
    double average;
    double sigma[3];
    vector<double> yValues; 

    for (int j = 0; j < 3; j++){
        yValues.clear(); // czyści wektor
        double sum = 0;
        double bigSum = 0;

        for (int i = 40 * j ; i < 40 * j + 40; i++){ // bierze pierwsze 40 punktów, potem kolejne 40
            double y = spectrum_ch->GetY()[i];
            yValues.push_back(y);
        }

        for (int i=0; i < yValues.size(); i++){
            sum += yValues[i];
        }
        average = sum / yValues.size();
        
        for(int i=0; i < yValues.size(); i++){
            bigSum += pow(yValues[i] - average, 2);
        }

        sigma[j] = sqrt(bigSum / (yValues.size() - 1));
    }

    double SmallestSigma = sigma[0];

    if (sigma[1] < SmallestSigma)
        SmallestSigma = sigma[1];

    if (sigma[2] < SmallestSigma)
        SmallestSigma = sigma[2];

    return SmallestSigma;
}

TGraph* EstimateBackground(vector<double> arguments_ch, vector<double> values_ch)
{
// Looking for local minima
  vector<int> minIndex;
  int minBinHalfRangeForMin = 5;
  for (int i=minBinHalfRangeForMin; i<arguments_ch.size()-minBinHalfRangeForMin; i++) {
    bool possibleMin = true;
    for (int j=1; j<=minBinHalfRangeForMin; j++) {
      if (values_ch.at(i) > values_ch.at(i-j) || values_ch.at(i) > values_ch.at(i+j))
        possibleMin = false;
    }
    if (possibleMin) {
      minIndex.push_back(i);
    }
  }
//
// Looking for standard deviation around the local minima
  vector<double> sdAroundMin;
  for (int i=0; i<minIndex.size(); i++) {
    double mean = 0;
    for (int j=-minBinHalfRangeForMin; j<=minBinHalfRangeForMin; j++) {
      mean += values_ch.at(minIndex.at(i)+j)/(2*minBinHalfRangeForMin+1);
    }
    double sd = 0;
    for (int j=-minBinHalfRangeForMin; j<=minBinHalfRangeForMin; j++) {
      sd += std::pow(values_ch.at(minIndex.at(i)+j) - mean, 2)/(2*minBinHalfRangeForMin+1)/(2*minBinHalfRangeForMin);
    }
    sdAroundMin.push_back(sqrt(sd));
  }
//
// Getting rid of sd too high and calculating mean SD
  vector<double> copyForMedian = sdAroundMin;
  sort(copyForMedian.begin(), copyForMedian.end());
  double median = copyForMedian.at((int)(0.5*copyForMedian.size()));
  double ceilingForSD = 2*median;
  double meanSD = 0;
  int noOfNumberToMean = 0;
  vector<double> minIndexCropped;
  for (unsigned i=0; i<sdAroundMin.size(); i++) {
    if (sdAroundMin.at(i) < ceilingForSD) {
      meanSD += sdAroundMin.at(i);
      noOfNumberToMean++;
      minIndexCropped.push_back(minIndex.at(i));
    }
  }
  meanSD /= noOfNumberToMean;
  cout << "Median and meanSD: " << median << " " << meanSD << endl;
  tempMeanSD = meanSD;
//
// Getting rid of minimal points that are coming from the peak defined as the difference between consecutive points greater than
// discardConstant * meanSD
  double discardConstant = 2.75;
  double discardValue = discardConstant*meanSD;
  vector<double> cleanerMinIndex;
  double sthLowerThanTheLowestPoint = -1000;
  double prevMin = sthLowerThanTheLowestPoint, prev2Min = sthLowerThanTheLowestPoint;
  for (unsigned i=0; i<minIndexCropped.size(); i++) {
    if (prevMin == sthLowerThanTheLowestPoint) {
      if (i == minIndexCropped.size()-1)
        break;
      if (fabs(values_ch.at(minIndexCropped.at(i)) - values_ch.at(minIndexCropped.at(i+1))) <= discardValue) {
// needs to be like that to cover the situation with random ultra low bins at the beginning
        prevMin = values_ch.at(minIndexCropped.at(i));
        cleanerMinIndex.push_back(minIndexCropped.at(i));
      }
    } else if (i<minIndexCropped.size()-1) {
      if (values_ch.at(minIndexCropped.at(i)) <= prevMin + discardValue &&
          values_ch.at(minIndexCropped.at(i)) <= values_ch.at(minIndexCropped.at(i+1)) + discardValue) {
        prev2Min = prevMin;
        prevMin = values_ch.at(minIndexCropped.at(i));
        cleanerMinIndex.push_back(minIndexCropped.at(i));
      } else if (prev2Min != sthLowerThanTheLowestPoint) {
// covering situation with some random ultra low bins
        if (values_ch.at(minIndexCropped.at(i)) <= prev2Min + discardValue &&
            values_ch.at(minIndexCropped.at(i)) <= values_ch.at(minIndexCropped.at(i+1)) + discardValue) {
          prev2Min = prevMin;
          prevMin = values_ch.at(minIndexCropped.at(i));
          cleanerMinIndex.push_back(minIndexCropped.at(i));
        }
      }
    } else {
      if (values_ch.at(minIndexCropped.at(i)) <= prevMin + discardValue) {
        cleanerMinIndex.push_back(minIndexCropped.at(i));
      }
    }
  }
  cout << "Index size: raw - cropped - cleaned" << minIndex.size() << " - " << minIndexCropped.size() << " - " << cleanerMinIndex.size() << endl;
//
// Drawing and returning background
  Double_t graphX[cleanerMinIndex.size()], graphY[cleanerMinIndex.size()];
  for (unsigned i=0; i<cleanerMinIndex.size(); i++) {
    graphX[i] = arguments_ch.at(cleanerMinIndex.at(i));
    graphY[i] = values_ch.at(cleanerMinIndex.at(i));
  }
  return new TGraph(cleanerMinIndex.size(), graphX, graphY);
}

TGraph* FindPeaks(vector<double> arguments_ch, vector<double> values_ch, TF1* bgFit)
{
  vector<double> peakArg, peakVal;
  double noiseRange = 12*tempMeanSD;
  unsigned counter = 0;
  unsigned minCounter = 2;
  int minCounterForPeak = 7;
  for (unsigned i=0; i<arguments_ch.size(); i++) {
    if (values_ch.at(i) > bgFit->Eval(arguments_ch.at(i)) + noiseRange) {
      counter++;
    } else {
      if (counter > minCounter) {
        int noOfPotentialPeaks = floor(counter/(minCounterForPeak+1)) + 1;
// minCounter + 1 to have a proper cut for multiple peaks
        int skipBetweenPeaks = floor(counter/(noOfPotentialPeaks+1));
        for (unsigned j=0; j<noOfPotentialPeaks; j++) {
          int currIt = i - (j+1)*skipBetweenPeaks;
          peakArg.push_back(arguments_ch.at(currIt));
          peakVal.push_back(values_ch.at(currIt) - (bgFit->Eval(arguments_ch.at(currIt)) + noiseRange));
        }
      }
      counter = 0;
    }
  }
  if (peakArg.size() == 0)
    return new TGraph();

  return new TGraph(peakArg.size(), peakArg.data(), peakVal.data());
}

Double_t peakFunc(Double_t *A, Double_t *P)
{
  Double_t result = 0;
  result += P[0] + P[1]*A[0] + P[2]*std::pow(A[0], 2) + P[3]*std::pow(A[0], 3) + P[4]*std::pow(A[0], 4);

  Int_t noOfRandomVoigt = P[5];
  Double_t gamma = P[6];
  Double_t sigma = P[7];
  for (unsigned i=0; i<noOfRandomVoigt; i++) {
    result += P[8+2*i]*VoigtForPeak(A[0], P[9+2*i], sigma, gamma);
  }

  return result;
}
