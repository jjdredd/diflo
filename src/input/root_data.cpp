#include "data.hpp"


//******************************************************************************
//*Tree    :out       : out                                                    *
//*Entries :   100000 : Total =      8982045806 bytes  File  Size = 4315392393 *
//*        :          : Tree compression factor =   2.08                       *
//******************************************************************************
//*Br    0 :npart     : npart/I                                                *
//*Entries :   100000 : Total  Size=     413580 bytes  File Size  =     195165 *
//*Baskets :      146 : Basket Size=      51200 bytes  Compression=   2.10     *
//*............................................................................*
//*Br    1 :x         : x[npart]/D                                             *
//*Entries :   100000 : Total  Size=  921051022 bytes  File Size  =  348712652 *
//*Baskets :     1435 : Basket Size=    4594176 bytes  Compression=   2.64     *
//*............................................................................*
//*Br    2 :y         : y[npart]/D                                             *
//*Entries :   100000 : Total  Size=  921051022 bytes  File Size  =  355313678 *
//*Baskets :     1435 : Basket Size=    4594176 bytes  Compression=   2.59     *
//*............................................................................*
//*Br    3 :z         : z[npart]/D                                             *
//*Entries :   100000 : Total  Size=  921051022 bytes  File Size  =  345459572 *
//*Baskets :     1435 : Basket Size=    4594176 bytes  Compression=   2.67     *
//*............................................................................*
//*Br    4 :t         : t[npart]/D                                             *
//*Entries :   100000 : Total  Size=  921051022 bytes  File Size  =  180003421 *
//*Baskets :     1435 : Basket Size=    4594176 bytes  Compression=   5.12     *
//*............................................................................*
//*Br    5 :px        : px[npart]/D                                            *
//*Entries :   100000 : Total  Size=  921052461 bytes  File Size  =  703361387 *
//*Baskets :     1435 : Basket Size=    4594688 bytes  Compression=   1.31     *
//*............................................................................*
//*Br    6 :py        : py[npart]/D                                            *
//*Entries :   100000 : Total  Size=  921052461 bytes  File Size  =  714814925 *
//*Baskets :     1435 : Basket Size=    4594688 bytes  Compression=   1.29     *
//*............................................................................*
//*Br    7 :pz        : pz[npart]/D                                            *
//*Entries :   100000 : Total  Size=  921052461 bytes  File Size  =  705973213 *
//*Baskets :     1435 : Basket Size=    4594688 bytes  Compression=   1.30     *
//*............................................................................*
//*Br    8 :E         : E[npart]/D                                             *
//*Entries :   100000 : Total  Size=  921051022 bytes  File Size  =  705997076 *
//*Baskets :     1435 : Basket Size=    4594176 bytes  Compression=   1.30     *
//*............................................................................*
//*Br    9 :id        : id[npart]/I                                            *
//*Entries :   100000 : Total  Size=  460737302 bytes  File Size  =   90029645 *
//*Baskets :      831 : Basket Size=    2296832 bytes  Compression=   5.12     *
//*............................................................................*
//*Br   10 :mid       : mid[npart]/I                                           *
//*Entries :   100000 : Total  Size=  460738137 bytes  File Size  =   77777613 *
//*Baskets :      831 : Basket Size=    2296832 bytes  Compression=   5.92     *
//*............................................................................*
//*Br   11 :ele       : ele[npart]/S                                           *
//*Entries :   100000 : Total  Size=  230581624 bytes  File Size  =   45319278 *
//*Baskets :      543 : Basket Size=    1148928 bytes  Compression=   5.09     *
//*............................................................................*
//*Br   12 :bar       : bar[npart]/S                                           *
//*Entries :   100000 : Total  Size=  230581624 bytes  File Size  =   22326899 *
//*Baskets :      543 : Basket Size=    1148928 bytes  Compression=  10.33     *
//*............................................................................*
//*Br   13 :str       : str[npart]/S                                           *
//*Entries :   100000 : Total  Size=  230581624 bytes  File Size  =   19962756 *
//*Baskets :      543 : Basket Size=    1148928 bytes  Compression=  11.55     *
//*............................................................................*
// None



ROOTData::ROOTData() {}

ROOTData::~ROOTData() {}

ROOTData::ROOTData(char * file_name) : tfin(file_name), cur_event(0) {
	tin = tfin.Get("out");
	nentries = tin->GetEntries();
	tin->Print();
}

bool ROOTData::FetchEvent(event &e) {

	unsigned npart;
	bool ret;

	e.particles.clear();

	tin->SetBranchStatus("*", 0);
	tin->SetBranchStatus("npart", 1);
	tin->SetBranchAddress("npart", &npart);
	if (tin->GetEntry(cur_event) <= 0) return false;

	e.reserve(npart);

	double *x = new double[npart],
		*y = new double[npart],
		*z = new double[npart],
		*px = new double[npart],
		*py = new double[npart],
		*pz = new double[npart],
		*E = new double[npart];
	int *id = new int[npart];

	tin->SetBranchStatus("*", 1);
	tin->SetBranchAddress("x", x);
	tin->SetBranchAddress("y", y);
	tin->SetBranchAddress("z", z);
	tin->SetBranchAddress("px", px);
	tin->SetBranchAddress("py", py);
	tin->SetBranchAddress("pz", pz);
	tin->SetBranchAddress("E", E);
	tin->SetBranchAddress("id", id);

	if (tin->GetEntry(cur_event) > 0) {
		// copy all that stuff
		particle p;
		for (unsigned i = 0; i < npart; i++) {
			p.type = id[i];
			p.Px = px[i];
			p.Py = py[i];
			p.Pz = pz[i];
			p.P0 = E[i];
			p.x = x[i];
			p.y = y[i];
			p.z = z[i];
			e.add_particle(p);
		}
		ret = true;
	} else ret = false;

	delete x;
	delete y;
	delete z;
	delete px;
	delete py;
	delete pz;
	delete E;
	delete id;

	return ret;
}
