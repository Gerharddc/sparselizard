// sparselizard - Copyright (C) see copyright file.
//
// See the LICENSE file for license information. Please report all
// bugs and problems to <alexandre.halbach at gmail.com>.


// This object manages the decomposition of the mesh into domains.

#ifndef DTRACKER_H
#define DTRACKER_H

#include <iostream>
#include <vector>
#include "rawmesh.h"
#include "slmpi.h"


class dtracker
{

    private:

        std::weak_ptr<rawmesh> myrawmesh;

        // Neighbours (without this rank and without duplicates, sorted):
        std::vector<int> myneighbours = {};
        // Direct access:
        std::vector<bool> myisneighbour = {};
        // No overlap interfaces:
        std::vector<int> mynooverlapinterfaces = {};


        // Discover up to 'numtrialelements' neighbours that share cell-1 dimension elements with this rank.
        // The barycenters of all element shared with the neighbours to discover must be provided as argument.
        // The number of interface elements on all ranks is returned (if all zero an empty vector is returned).
        std::vector<int> discoversomeneighbours(int numtrialelements, std::vector<double>& interfaceelembarys, std::vector<int>& neighboursfound);

        // Upon return 'inneighbours[i]' is the number of the neighbour touching the ith interface element (-1 if no neighbour touching).
        // The neighbours provided must be unique. The number of interface elements for each rank must be provided in 'allnumelementsininterface'.
        void discoverinterfaces(std::vector<int> neighbours, std::vector<double>& interfaceelembarys, std::vector<int>& allnumelementsininterface, std::vector<int>& inneighbour);

        // Find new interfaces and populate the output accordingly. Return false if no more interface can be found on any rank.
        bool discovercrossinterfaces(std::vector<int>& interfacenodelist, std::vector<int>& interfaceedgelist, std::vector<std::vector<bool>>& isnodeinneighbours, std::vector<std::vector<bool>>& isedgeinneighbours);

    public:

        dtracker(std::shared_ptr<rawmesh> rm);

        std::shared_ptr<rawmesh> getrawmesh(void);
        
        // Set manually or discover automatically the no-overlap connectivity of this rank:
        void setconnectivity(std::vector<int> neighbours, std::vector<int> nooverlapinterfaces);
        void discoverconnectivity(int nooverlapinterface, int numtrialelements = 10, int verbosity = 0);

        int countneighbours(void);
        int getneighbour(int neighbourindex);

        bool isneighbour(int neighbour);
        int getnooverlapinterface(int neighbour);

        // Print connectivity information:
        void print(void);
        // Write no-overlap domain interfaces:
        void writeinterfaces(std::string filename);

};

#endif
