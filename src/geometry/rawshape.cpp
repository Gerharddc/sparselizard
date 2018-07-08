#include "rawshape.h"
#include "geotools.h"


void rawshape::deform(expression xdeform, expression ydeform, expression zdeform, bool recursively)
{
	std::vector<double>* mycoords = getcoords();

	int numnodes = mycoords->size()/3;

	std::vector<double> xcoords(numnodes);
	std::vector<double> ycoords(numnodes);
	std::vector<double> zcoords(numnodes);

	for (int i = 0; i < numnodes; i++)
	{
		xcoords[i] = mycoords->at(3*i+0);
		ycoords[i] = mycoords->at(3*i+1);
		zcoords[i] = mycoords->at(3*i+2);
	}

	std::vector<double> xdeformvec = xdeform.evaluate(xcoords,ycoords,zcoords);
	std::vector<double> ydeformvec = ydeform.evaluate(xcoords,ycoords,zcoords);
	std::vector<double> zdeformvec = zdeform.evaluate(xcoords,ycoords,zcoords);

	for (int i = 0; i < numnodes; i++)
	{
		mycoords->at(3*i+0) += xdeformvec[i];
		mycoords->at(3*i+1) += ydeformvec[i];
		mycoords->at(3*i+2) += zdeformvec[i];
	}

	
	// Also deform (only once) every sub shape:
	if (recursively)
	{
		std::vector<rawshape*> subshapes = geotools::unique(geotools::getpointers(getsubshapesrecursively()));
		for (int i = 0; i < subshapes.size(); i++)
			subshapes[i]->deform(xdeform, ydeform, zdeform, false);
	}
}

void rawshape::shift(double shiftx, double shifty, double shiftz, bool recursively)
{
	std::vector<double>* mycoords = getcoords();

	int numnodes = mycoords->size()/3;

	for (int i = 0; i < numnodes; i++)
	{
		mycoords->at(3*i+0) += shiftx;
		mycoords->at(3*i+1) += shifty;
		mycoords->at(3*i+2) += shiftz;
	}


	// Also shift (only once) every sub shape:
	if (recursively)
	{
		std::vector<rawshape*> subshapes = geotools::unique(geotools::getpointers(getsubshapesrecursively()));
		for (int i = 0; i < subshapes.size(); i++)
			subshapes[i]->shift(shiftx, shifty, shiftz, false);
	}
}

void rawshape::rotate(double alphax, double alphay, double alphaz, bool recursively)
{
	// Get the coordinates of the shape to rotate:
	std::vector<double>* mycoords = getcoords();

    geotools::rotate(alphax, alphay, alphaz, mycoords);


	// Also rotate (only once) every sub shape:
	if (recursively)
	{
		std::vector<rawshape*> subshapes = geotools::unique(geotools::getpointers(getsubshapesrecursively()));
		for (int i = 0; i < subshapes.size(); i++)
			subshapes[i]->rotate(alphax, alphay, alphaz, false);
	}
}


std::shared_ptr<rawshape> rawshape::extrude(int physreg, double height, int numlayers)
{
	std::cout << "Error in 'rawshape' object: 'extrude' has not been defined for this shape" << std::endl;
	abort(); 
}

std::shared_ptr<rawshape> rawshape::duplicate(void)
{
	std::cout << "Error in 'rawshape' object: 'duplicate' has not been defined for this shape" << std::endl;
	abort(); 
}

void rawshape::flip(void)
{
	std::cout << "Error in 'rawshape' object: 'flip' has not been defined for this shape" << std::endl;
	abort(); 
}

void rawshape::setphysicalregion(int physreg) 
{
	std::cout << "Error in 'rawshape' object: 'setphysicalregion' has not been defined for this shape" << std::endl;
	abort(); 
}

int rawshape::getdimension(void)
{
	std::cout << "Error in 'rawshape' object: 'getdimension' has not been defined for this shape" << std::endl;
	abort(); 
}

std::string rawshape::getname(void)
{
	std::cout << "Error in 'rawshape' object: 'getname' has not been defined for this shape" << std::endl;
	abort(); 
}

std::vector<std::shared_ptr<rawshape>> rawshape::getsons(void)
{
	std::cout << "Error in 'rawshape' object: 'getsons' has not been defined for this shape" << std::endl;
	abort(); 
}

std::vector<std::shared_ptr<rawshape>> rawshape::getsubshapes(void)
{
	std::cout << "Error in 'rawshape' object: 'getsubshapes' has not been defined for this shape" << std::endl;
	abort(); 
}

std::vector<std::shared_ptr<rawshape>> rawshape::getsubshapesrecursively(void)
{
	std::vector<std::shared_ptr<rawshape>> subshapes = getsubshapes();

	if (subshapes.size() == 0)
		return subshapes;

	std::vector<std::shared_ptr<rawshape>> output = {};

	for (int i = 0; i < subshapes.size(); i++)
	{
		output.push_back(subshapes[i]);

		std::vector<std::shared_ptr<rawshape>> toappend = subshapes[i]->getsubshapesrecursively();
		for (int j = 0; j < toappend.size(); j++)
			output.push_back(toappend[j]);
	}

	return output;
}

int rawshape::getphysicalregion(void)
{
	std::cout << "Error in 'rawshape' object: 'getphysicalregion' has not been defined for this shape" << std::endl;
	abort(); 
}

std::vector<double>* rawshape::getcoords(void)
{
	std::cout << "Error in 'rawshape' object: 'getcoords' has not been defined for this shape" << std::endl;
	abort(); 
}

std::vector<std::vector<int>>* rawshape::getelems(void)
{
	std::cout << "Error in 'rawshape' object: 'getelems' has not been defined for this shape" << std::endl;
	abort(); 
}

std::shared_ptr<rawshape> rawshape::getpointer(void)
{
	std::cout << "Error in 'rawshape' object: 'getpointer' has not been defined for this shape" << std::endl;
	abort();
}

void rawshape::mesh(void)
{
	std::cout << "Error in 'rawshape' object: 'mesh' has not been defined for this shape" << std::endl;
	abort();
}



