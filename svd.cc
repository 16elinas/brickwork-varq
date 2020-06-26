//
// 2019 Many Electron Collaboration Summer School
// ITensor Tutorial
//
#include "itensor/all.h"

using namespace itensor;

int main()
    {
    //
    // SVD of CZ
    //

    auto i1 = Index(2,"i1");
    auto i2 = Index(2,"i2");
    auto ii1 = Index(2, "ii1");
    auto ii2 = Index(2, "ii2");

    auto cz = ITensor(i1, i2, ii1, ii2);
    Print(cz);

    auto [Ci, ci] = combiner(i1, i2);
    auto [Cii, cii] = combiner(ii1, ii2);

    auto Ccz = Ci * cz;
    Print(dim(ci));

    Ccz.set(ci=1, cii=1, 1);
    Ccz.set(ci=2, cii=2, 1);
    Ccz.set(ci=3, cii=3, 1);
    Ccz.set(ci=4, cii=4, -1);
    PrintData(Ccz);

    return 0;
    }
