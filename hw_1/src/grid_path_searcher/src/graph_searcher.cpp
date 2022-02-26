#include <graph_searcher.h>

using namespace std;
using namespace Eigen;

void gridPathFinder::initGridMap(double _resolution, Vector3d global_xyz_l, Vector3d global_xyz_u, int max_x_id, int max_y_id, int max_z_id)
{   
    gl_xl = global_xyz_l(0);
    gl_yl = global_xyz_l(1);
    gl_zl = global_xyz_l(2);

    gl_xu = global_xyz_u(0);
    gl_yu = global_xyz_u(1);
    gl_zu = global_xyz_u(2);
    
    GLX_SIZE = max_x_id;
    GLY_SIZE = max_y_id;
    GLZ_SIZE = max_z_id;
    GLYZ_SIZE  = GLY_SIZE * GLZ_SIZE;
    GLXYZ_SIZE = GLX_SIZE * GLYZ_SIZE;

    resolution = _resolution;
    inv_resolution = 1.0 / _resolution;    

    data = new uint8_t[GLXYZ_SIZE];
    memset(data, 0, GLXYZ_SIZE * sizeof(uint8_t));
    
    GridNodeMap = new GridNodePtr ** [GLX_SIZE];
    for(int i = 0; i < GLX_SIZE; i++){
        GridNodeMap[i] = new GridNodePtr * [GLY_SIZE];
        for(int j = 0; j < GLY_SIZE; j++){
            GridNodeMap[i][j] = new GridNodePtr [GLZ_SIZE];
            for( int k = 0; k < GLZ_SIZE;k++){
                Vector3i tmpIdx(i,j,k);
                Vector3d pos = gridIndex2coord(tmpIdx);
                GridNodeMap[i][j][k] = new GridNode(tmpIdx, pos);
            }
        }
    }
}

void gridPathFinder::resetGrid(GridNodePtr ptr)
{
    ptr->id = 0;
    ptr->cameFrom = NULL;
    ptr->gScore = inf;
    ptr->fScore = inf;
}

void gridPathFinder::resetUsedGrids()
{   
    //ROS_WARN("expandedNodes size : %d", expandedNodes.size()); 
    for(auto tmpPtr:expandedNodes)
        resetGrid(tmpPtr);

    GridNodePtr tmpPtr = NULL;
    for(auto ptr:openSet){   
        tmpPtr = ptr.second;
        resetGrid(tmpPtr);
    }

    for(auto ptr:gridPath)
        ptr->is_path = false;

    expandedNodes.clear();
}

void gridPathFinder::setObs(const double coord_x, const double coord_y, const double coord_z)
{   
    if( coord_x < gl_xl  || coord_y < gl_yl  || coord_z <  gl_zl || 
        coord_x >= gl_xu || coord_y >= gl_yu || coord_z >= gl_zu )
        return;

    int idx_x = static_cast<int>( (coord_x - gl_xl) * inv_resolution);
    int idx_y = static_cast<int>( (coord_y - gl_yl) * inv_resolution);
    int idx_z = static_cast<int>( (coord_z - gl_zl) * inv_resolution);      

    data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] = 1;
}

double gridPathFinder::getDiagHeu(GridNodePtr node1, GridNodePtr node2)
{   
    double dx = abs(node1->index(0) - node2->index(0));
    double dy = abs(node1->index(1) - node2->index(1));
    double dz = abs(node1->index(2) - node2->index(2));

    double h = 0.0;
    int diag = min(min(dx, dy), dz);
    dx -= diag;
    dy -= diag;
    dz -= diag;

    if (dx == 0) {
        h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dy, dz) + 1.0 * abs(dy - dz);
    }
    if (dy == 0) {
        h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dx, dz) + 1.0 * abs(dx - dz);
    }
    if (dz == 0) {
        h = 1.0 * sqrt(3.0) * diag + sqrt(2.0) * min(dx, dy) + 1.0 * abs(dx - dy);
    }
    return h;
}

double gridPathFinder::getManhHeu(GridNodePtr node1, GridNodePtr node2)
{   
    double dx = abs(node1->index(0) - node2->index(0));
    double dy = abs(node1->index(1) - node2->index(1));
    double dz = abs(node1->index(2) - node2->index(2));

    return dx + dy + dz;
}

double gridPathFinder::getEuclHeu(GridNodePtr node1, GridNodePtr node2)
{   
    return (node2->index - node1->index).norm();
}

double gridPathFinder::getHeu(GridNodePtr node1, GridNodePtr node2)
{
    return tie_breaker * getDiagHeu(node1, node2);
    //return getEuclHeu(node1, node2);
}

vector<GridNodePtr> gridPathFinder::retrievePath(GridNodePtr current)
{   
    vector<GridNodePtr> path;
    path.push_back(current);

    while(current->cameFrom != NULL){   
        current->is_path = true;
        current = current -> cameFrom;
        path.push_back(current);
    }

    current->is_path = true;

    return path;
}

vector<Vector3d> gridPathFinder::getVisitedNodes()
{   
    vector<Vector3d> visited_nodes;
    for(int i = 0; i < GLX_SIZE; i++)
        for(int j = 0; j < GLY_SIZE; j++)
            for(int k = 0; k < GLZ_SIZE; k++){   
                if(GridNodeMap[i][j][k]->id != 0)
                //if(GridNodeMap[i][j][k]->id == -1)
                    visited_nodes.push_back(GridNodeMap[i][j][k]->coord);
            }

    ROS_WARN("visited_nodes size : %d", visited_nodes.size());
    return visited_nodes;
}

vector<Vector3d> gridPathFinder::getCloseNodes()
{   
    vector<Vector3d> vec;
    for(auto tmpPtr:expandedNodes){   
        if( !tmpPtr->is_path )
            vec.push_back(tmpPtr->coord);
    }

    return vec;
}

vector<Vector3d> gridPathFinder::getPath() 
{   
    vector<Vector3d> path;

    for(auto ptr: gridPath)
        path.push_back(ptr->coord);

    reverse(path.begin(), path.end());
    return path;
}

inline Vector3d gridPathFinder::gridIndex2coord(const Vector3i & index) const
{
    Vector3d pt;

    pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
    pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
    pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;

    return pt;
}

inline Vector3i gridPathFinder::coord2gridIndex(const Vector3d & pt) const
{
    Vector3i idx;
    idx <<  min( max( int( (pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1),
            min( max( int( (pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1),
            min( max( int( (pt(2) - gl_zl) * inv_resolution), 0), GLZ_SIZE - 1);                  
  
    return idx;
}

Eigen::Vector3d gridPathFinder::coordRounding(const Eigen::Vector3d & coord) const
{
    return gridIndex2coord(coord2gridIndex(coord));
}

inline bool gridPathFinder::isOccupied(const Eigen::Vector3i & index) const
{
    return isOccupied(index(0), index(1), index(2));
}

inline bool gridPathFinder::isFree(const Eigen::Vector3i & index) const
{
    return isFree(index(0), index(1), index(2));
}

inline bool gridPathFinder::isOccupied(const int & idx_x, const int & idx_y, const int & idx_z) const 
{
    return  (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE && 
            (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] == 1));
}

inline bool gridPathFinder::isFree(const int & idx_x, const int & idx_y, const int & idx_z) const 
{
    return (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE && 
           (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] < 1));
}

inline void gridPathFinder::getJpsSucc(GridNodePtr currentPtr, vector<GridNodePtr> & neighborPtrSets, vector<double> & edgeCostSets, int num_iter)
{   
    neighborPtrSets.clear();
    edgeCostSets.clear();
    const int norm1 = abs(currentPtr->dir(0)) + abs(currentPtr->dir(1)) + abs(currentPtr->dir(2));

    int num_neib  = jn3d->nsz[norm1][0];
    int num_fneib = jn3d->nsz[norm1][1];
    int id = (currentPtr->dir(0) + 1) + 3 * (currentPtr->dir(1) + 1) + 9 * (currentPtr->dir(2) + 1);

    for( int dev = 0; dev < num_neib + num_fneib; ++dev) {
        Vector3i neighborIdx;
        Vector3i expandDir;

        if( dev < num_neib ) {
            expandDir(0) = jn3d->ns[id][0][dev];
            expandDir(1) = jn3d->ns[id][1][dev];
            expandDir(2) = jn3d->ns[id][2][dev];
            
            if( !jump(currentPtr->index, expandDir, neighborIdx) )  
                continue;
        }
        else {
            int nx = currentPtr->index(0) + jn3d->f1[id][0][dev - num_neib];
            int ny = currentPtr->index(1) + jn3d->f1[id][1][dev - num_neib];
            int nz = currentPtr->index(2) + jn3d->f1[id][2][dev - num_neib];
            
            if( isOccupied(nx, ny, nz) ) {
                expandDir(0) = jn3d->f2[id][0][dev - num_neib];
                expandDir(1) = jn3d->f2[id][1][dev - num_neib];
                expandDir(2) = jn3d->f2[id][2][dev - num_neib];
                
                if( !jump(currentPtr->index, expandDir, neighborIdx) ) 
                    continue;
            }
            else
                continue;
        }

        if( num_iter == 1 )
            debugNodes.push_back(gridIndex2coord(neighborIdx));

        GridNodePtr nodePtr = GridNodeMap[neighborIdx(0)][neighborIdx(1)][neighborIdx(2)];
        nodePtr->dir = expandDir;
        
        neighborPtrSets.push_back(nodePtr);
        edgeCostSets.push_back(
            sqrt(
            (neighborIdx(0) - currentPtr->index(0)) * (neighborIdx(0) - currentPtr->index(0)) +
            (neighborIdx(1) - currentPtr->index(1)) * (neighborIdx(1) - currentPtr->index(1)) +
            (neighborIdx(2) - currentPtr->index(2)) * (neighborIdx(2) - currentPtr->index(2))   ) 
            );
    }
}

inline void gridPathFinder::getSucc(GridNodePtr currentPtr, vector<GridNodePtr> & neighborPtrSets, vector<double> & edgeCostSets)
{   
    neighborPtrSets.clear();
    edgeCostSets.clear();
    Vector3i neighborIdx;
    for(int dx = -1; dx < 2; dx++){
        for(int dy = -1; dy < 2; dy++){
            for(int dz = -1; dz < 2; dz++){
                
                if( dx == 0 && dy == 0 && dz ==0 )
                    continue; 

                neighborIdx(0) = (currentPtr -> index)(0) + dx;
                neighborIdx(1) = (currentPtr -> index)(1) + dy;
                neighborIdx(2) = (currentPtr -> index)(2) + dz;

                if(    neighborIdx(0) < 0 || neighborIdx(0) >= GLX_SIZE
                    || neighborIdx(1) < 0 || neighborIdx(1) >= GLY_SIZE
                    || neighborIdx(2) < 0 || neighborIdx(2) >= GLZ_SIZE){
                    continue;
                }

                neighborPtrSets.push_back(GridNodeMap[neighborIdx(0)][neighborIdx(1)][neighborIdx(2)]);
                edgeCostSets.   push_back(sqrt(dx * dx + dy * dy + dz * dz));
            }
        }
    }
}

void gridPathFinder::graphSearch(Vector3d start_pt, Vector3d end_pt, bool use_jps)
{   
    ros::Time time_1 = ros::Time::now();    
    debugNodes.clear();

    Vector3i start_idx = coord2gridIndex(start_pt);
    Vector3i end_idx   = coord2gridIndex(end_pt);

    goalIdx = end_idx;

    start_pt = gridIndex2coord(start_idx);
    end_pt   = gridIndex2coord(end_idx);

    GridNodePtr startPtr = new GridNode(start_idx, start_pt);
    GridNodePtr endPtr   = new GridNode(end_idx,   end_pt);

    openSet.clear();

    GridNodePtr neighborPtr = NULL;
    GridNodePtr currentPtr  = NULL;

    startPtr -> gScore = 0;
    startPtr -> fScore = getHeu(startPtr, endPtr);
    startPtr -> id = 1; //put start node in open set
    startPtr -> coord = start_pt;
    openSet.insert( make_pair(startPtr -> fScore, startPtr) ); //put start in open set

    double tentative_gScore;

    int num_iter = 0;
    vector<GridNodePtr> neighborPtrSets;
    vector<double> edgeCostSets;

    // we only cover 3d case in this project.
    while ( !openSet.empty() ){   
        num_iter ++;
        currentPtr = openSet.begin() -> second;

        if( currentPtr->index == goalIdx ){
            ros::Time time_2 = ros::Time::now();

            if( use_jps )
                ROS_WARN("[JPS]{sucess} Time in JPS is %f ms, path cost if %f m", (time_2 - time_1).toSec() * 1000.0, currentPtr->gScore * resolution );
            else
                ROS_WARN("[A*]{sucess}  Time in A*  is %f ms, path cost if %f m", (time_2 - time_1).toSec() * 1000.0, currentPtr->gScore * resolution );
            
            gridPath = retrievePath(currentPtr);
            return;
        }         
        openSet.erase(openSet.begin());
        currentPtr -> id = -1; //move current node from open set to closed set.
        expandedNodes.push_back(currentPtr);
        
        if(!use_jps)
            getSucc(currentPtr, neighborPtrSets, edgeCostSets);
        else
            getJpsSucc(currentPtr, neighborPtrSets, edgeCostSets, num_iter);

        for(int i = 0; i < (int)neighborPtrSets.size(); i++){
            neighborPtr = neighborPtrSets[i];
            if( isOccupied(neighborPtr->index) || neighborPtr -> id == -1 )
                continue;

            double edge_cost = edgeCostSets[i];            
            tentative_gScore = currentPtr -> gScore + edge_cost; 

            if(neighborPtr -> id != 1){ //discover a new node
                neighborPtr -> id        = 1;
                neighborPtr -> cameFrom  = currentPtr;
                neighborPtr -> gScore    = tentative_gScore;
                neighborPtr -> fScore    = neighborPtr -> gScore + getHeu(neighborPtr, endPtr); 
                neighborPtr -> nodeMapIt = openSet.insert( make_pair(neighborPtr->fScore, neighborPtr) ); //put neighbor in open set and record it.
                continue;
            }
            else if(tentative_gScore <= neighborPtr-> gScore){ //in open set and need update
                neighborPtr -> cameFrom = currentPtr;
                neighborPtr -> gScore = tentative_gScore;
                neighborPtr -> fScore = tentative_gScore + getHeu(neighborPtr, endPtr); 
                openSet.erase(neighborPtr -> nodeMapIt);
                neighborPtr -> nodeMapIt = openSet.insert( make_pair(neighborPtr->fScore, neighborPtr) ); //put neighbor in open set and record it.

                // if change its parents, update the expanding direction
                for(int i = 0; i < 3; i++){
                    neighborPtr->dir(i) = neighborPtr->index(i) - currentPtr->index(i);
                    if( neighborPtr->dir(i) != 0)
                        neighborPtr->dir(i) /= abs( neighborPtr->dir(i) );
                }
            }      
        }
    }

    ros::Time time_2 = ros::Time::now();

    if((time_2 - time_1).toSec() > 0.1)
        ROS_WARN("Time consume in A star path finding is %f", (time_2 - time_1).toSec() );
}

bool gridPathFinder::jump(const Vector3i & curIdx, const Vector3i & expDir, Vector3i & neiIdx)
{
    neiIdx = curIdx + expDir;

    if( !isFree(neiIdx) )
        return false;

    if( neiIdx == goalIdx )
        return true;

    if( hasForced(neiIdx, expDir) )
        return true;

    const int id = (expDir(0) + 1) + 3 * (expDir(1) + 1) + 9 * (expDir(2) + 1);
    const int norm1 = abs(expDir(0)) + abs(expDir(1)) + abs(expDir(2));
    int num_neib = jn3d->nsz[norm1][0];

    for( int k = 0; k < num_neib - 1; ++k ){
        Vector3i newNeiIdx;
        Vector3i newDir(jn3d->ns[id][0][k], jn3d->ns[id][1][k], jn3d->ns[id][2][k]);
        if( jump(neiIdx, newDir, newNeiIdx) ) 
            return true;
    }

    return jump(neiIdx, expDir, neiIdx);
}

inline bool gridPathFinder::hasForced(const Vector3i & idx, const Vector3i & dir)
{
    int norm1 = abs(dir(0)) + abs(dir(1)) + abs(dir(2));
    int id    = (dir(0) + 1) + 3 * (dir(1) + 1) + 9 * (dir(2) + 1);

    switch(norm1){
        case 1:
            // 1-d move, check 8 neighbors
            for( int fn = 0; fn < 8; ++fn ){
                int nx = idx(0) + jn3d->f1[id][0][fn];
                int ny = idx(1) + jn3d->f1[id][1][fn];
                int nz = idx(2) + jn3d->f1[id][2][fn];
                if( isOccupied(nx, ny, nz) )
                    return true;
            }
            return false;

        case 2:
            // 2-d move, check 8 neighbors
            for( int fn = 0; fn < 8; ++fn ){
                int nx = idx(0) + jn3d->f1[id][0][fn];
                int ny = idx(1) + jn3d->f1[id][1][fn];
                int nz = idx(2) + jn3d->f1[id][2][fn];
                if( isOccupied(nx, ny, nz) )
                    return true;
            }
            return false;

        case 3:
            // 3-d move, check 6 neighbors
            for( int fn = 0; fn < 6; ++fn ){
                int nx = idx(0) + jn3d->f1[id][0][fn];
                int ny = idx(1) + jn3d->f1[id][1][fn];
                int nz = idx(2) + jn3d->f1[id][2][fn];
                if( isOccupied(nx, ny, nz) )
                    return true;
            }
            return false;

        default:
            return false;
    }
}

constexpr int JPS3DNeib::nsz[4][2];
JPS3DNeib::JPS3DNeib() 
{
    int id = 0;
    for(int dz = -1; dz <= 1; ++ dz) {
        for(int dy = -1; dy <= 1; ++ dy) {
            for(int dx = -1; dx <= 1; ++ dx) {
                int norm1 = abs(dx) + abs(dy) + abs(dz);
            
                for(int dev = 0; dev < nsz[norm1][0]; ++ dev)
                    Neib(dx,dy,dz,norm1,dev, ns[id][0][dev], ns[id][1][dev], ns[id][2][dev]);
            
                for(int dev = 0; dev < nsz[norm1][1]; ++ dev){
                    FNeib(dx,dy,dz,norm1,dev,
                    f1[id][0][dev],f1[id][1][dev], f1[id][2][dev],
                    f2[id][0][dev],f2[id][1][dev], f2[id][2][dev]);
                }
                
                id ++;
            }
        }
    }
}


void JPS3DNeib::Neib(int dx, int dy, int dz, int norm1, int dev,
    int& tx, int& ty, int& tz)
{
    switch(norm1){
        case 0:
            switch(dev){
                case 0:  tx=1;  ty=0;  tz=0;  return;
                case 1:  tx=-1; ty=0;  tz=0;  return;
                case 2:  tx=0;  ty=1;  tz=0;  return;
                case 3:  tx=1;  ty=1;  tz=0;  return;
                case 4:  tx=-1; ty=1;  tz=0;  return;
                case 5:  tx=0;  ty=-1; tz=0;  return;
                case 6:  tx=1;  ty=-1; tz=0;  return;
                case 7:  tx=-1; ty=-1; tz=0;  return;
                case 8:  tx=0;  ty=0;  tz=1;  return;
                case 9:  tx=1;  ty=0;  tz=1;  return;
                case 10: tx=-1; ty=0;  tz=1;  return;
                case 11: tx=0;  ty=1;  tz=1;  return;
                case 12: tx=1;  ty=1;  tz=1;  return;
                case 13: tx=-1; ty=1;  tz=1;  return;
                case 14: tx=0;  ty=-1; tz=1;  return;
                case 15: tx=1;  ty=-1; tz=1;  return;
                case 16: tx=-1; ty=-1; tz=1;  return;
                case 17: tx=0;  ty=0;  tz=-1; return;
                case 18: tx=1;  ty=0;  tz=-1; return;
                case 19: tx=-1; ty=0;  tz=-1; return;
                case 20: tx=0;  ty=1;  tz=-1; return;
                case 21: tx=1;  ty=1;  tz=-1; return;
                case 22: tx=-1; ty=1;  tz=-1; return;
                case 23: tx=0;  ty=-1; tz=-1; return;
                case 24: tx=1;  ty=-1; tz=-1; return;
                case 25: tx=-1; ty=-1; tz=-1; return;
            }
        case 1:
            tx = dx; ty = dy; tz = dz; return;
        case 2:
            switch(dev){
                case 0:
                    if(dz == 0){
                        tx = 0; ty = dy; tz = 0; return;
                    }else{
                        tx = 0; ty = 0; tz = dz; return;
                    }
                case 1:
                    if(dx == 0){
                        tx = 0; ty = dy; tz = 0; return;
                    }else{
                        tx = dx; ty = 0; tz = 0; return;
                    }
                case 2:
                    tx = dx; ty = dy; tz = dz; return;
            }
        case 3:
            switch(dev){
                case 0: tx = dx; ty =  0; tz =  0; return;
                case 1: tx =  0; ty = dy; tz =  0; return;
                case 2: tx =  0; ty =  0; tz = dz; return;
                case 3: tx = dx; ty = dy; tz =  0; return;
                case 4: tx = dx; ty =  0; tz = dz; return;
                case 5: tx =  0; ty = dy; tz = dz; return;
                case 6: tx = dx; ty = dy; tz = dz; return;
            }
    }
}

void JPS3DNeib::FNeib( int dx, int dy, int dz, int norm1, int dev,
                          int& fx, int& fy, int& fz,
                          int& nx, int& ny, int& nz)
{
    switch(norm1){
        case 1:
            switch(dev){
                case 0: fx= 0; fy= 1; fz = 0; break;
                case 1: fx= 0; fy=-1; fz = 0; break;
                case 2: fx= 1; fy= 0; fz = 0; break;
                case 3: fx= 1; fy= 1; fz = 0; break;
                case 4: fx= 1; fy=-1; fz = 0; break;
                case 5: fx=-1; fy= 0; fz = 0; break;
                case 6: fx=-1; fy= 1; fz = 0; break;
                case 7: fx=-1; fy=-1; fz = 0; break;
            }
            nx = fx; ny = fy; nz = dz;
            // switch order if different direction
            if(dx != 0){
                fz = fx; fx = 0;
                nz = fz; nx = dx;
            }

            if(dy != 0){
                fz = fy; fy = 0;
                nz = fz; ny = dy;
            }
            return;
        case 2:
            if(dx == 0){
                switch(dev){
                    case 0:
                        fx = 0; fy = 0; fz = -dz;
                        nx = 0; ny = dy; nz = -dz;
                        return;
                    case 1:
                        fx = 0; fy = -dy; fz = 0;
                        nx = 0; ny = -dy; nz = dz;
                        return;
                    case 2:
                        fx = 1; fy = 0; fz = 0;
                        nx = 1; ny = dy; nz = dz;
                        return;
                    case 3:
                        fx = -1; fy = 0; fz = 0;
                        nx = -1; ny = dy; nz = dz;
                        return;
                    case 4:
                        fx = 1; fy = 0; fz = -dz;
                        nx = 1; ny = dy; nz = -dz;
                        return;
                    case 5:
                        fx = 1; fy = -dy; fz = 0;
                        nx = 1; ny = -dy; nz = dz;
                        return;
                    case 6:
                        fx = -1; fy = 0; fz = -dz;
                        nx = -1; ny = dy; nz = -dz;
                        return;
                    case 7:
                        fx = -1; fy = -dy; fz = 0;
                        nx = -1; ny = -dy; nz = dz;
                        return;
                    // Extras
                    case 8:
                        fx = 1; fy = 0; fz = 0;
                        nx = 1; ny = dy; nz = 0;
                        return;
                    case 9:
                        fx = 1; fy = 0; fz = 0;
                        nx = 1; ny = 0; nz = dz;
                        return;
                    case 10:
                        fx = -1; fy = 0; fz = 0;
                        nx = -1; ny = dy; nz = 0;
                        return;
                    case 11:
                        fx = -1; fy = 0; fz = 0;
                        nx = -1; ny = 0; nz = dz;
                        return;
                }
            }
            else if(dy == 0){
                switch(dev){
                    case 0:
                        fx = 0; fy = 0; fz = -dz;
                        nx = dx; ny = 0; nz = -dz;
                        return;
                    case 1:
                        fx = -dx; fy = 0; fz = 0;
                        nx = -dx; ny = 0; nz = dz;
                        return;
                    case 2:
                        fx = 0; fy = 1; fz = 0;
                        nx = dx; ny = 1; nz = dz;
                        return;
                    case 3:
                        fx = 0; fy = -1; fz = 0;
                        nx = dx; ny = -1;nz = dz;
                        return;
                    case 4:
                        fx = 0; fy = 1; fz = -dz;
                        nx = dx; ny = 1; nz = -dz;
                        return;
                    case 5:
                        fx = -dx; fy = 1; fz = 0;
                        nx = -dx; ny = 1; nz = dz;
                        return;
                    case 6:
                        fx = 0; fy = -1; fz = -dz;
                        nx = dx; ny = -1; nz = -dz;
                        return;
                    case 7:
                        fx = -dx; fy = -1; fz = 0;
                        nx = -dx; ny = -1; nz = dz;
                        return;
                    // Extras
                    case 8:
                        fx = 0; fy = 1; fz = 0;
                        nx = dx; ny = 1; nz = 0;
                        return;
                    case 9:
                        fx = 0; fy = 1; fz = 0;
                        nx = 0; ny = 1; nz = dz;
                        return;
                    case 10:
                        fx = 0; fy = -1; fz = 0;
                        nx = dx; ny = -1; nz = 0;
                        return;
                    case 11:
                        fx = 0; fy = -1; fz = 0;
                        nx = 0; ny = -1; nz = dz;
                        return;
                }
            }
            else{// dz==0
                switch(dev){
                    case 0:
                        fx = 0; fy = -dy; fz = 0;
                        nx = dx; ny = -dy; nz = 0;
                        return;
                    case 1:
                        fx = -dx; fy = 0; fz = 0;
                        nx = -dx; ny = dy; nz = 0;
                        return;
                    case 2:
                        fx =  0; fy = 0; fz = 1;
                        nx = dx; ny = dy; nz = 1;
                        return;
                    case 3:
                        fx =  0; fy = 0; fz = -1;
                        nx = dx; ny = dy; nz = -1;
                        return;
                    case 4:
                        fx = 0; fy = -dy; fz = 1;
                        nx = dx; ny = -dy; nz = 1;
                        return;
                    case 5:
                        fx = -dx; fy = 0; fz = 1;
                        nx = -dx; ny = dy; nz = 1;
                        return;
                    case 6:
                        fx = 0; fy = -dy; fz = -1;
                        nx = dx; ny = -dy; nz = -1;
                        return;
                    case 7:
                        fx = -dx; fy = 0; fz = -1;
                        nx = -dx; ny = dy; nz = -1;
                        return;
                    // Extras
                    case 8:
                        fx =  0; fy = 0; fz = 1;
                        nx = dx; ny = 0; nz = 1;
                        return;
                    case 9:
                        fx = 0; fy = 0; fz = 1;
                        nx = 0; ny = dy; nz = 1;
                        return;
                    case 10:
                        fx =  0; fy = 0; fz = -1;
                        nx = dx; ny = 0; nz = -1;
                        return;
                    case 11:
                        fx = 0; fy = 0; fz = -1;
                        nx = 0; ny = dy; nz = -1;
                        return;
                }
            }
        case 3:
            switch(dev){
                case 0:
                    fx = -dx; fy = 0; fz = 0;
                    nx = -dx; ny = dy; nz = dz;
                    return;
                case 1:
                    fx = 0; fy = -dy; fz = 0;
                    nx = dx; ny = -dy; nz = dz;
                    return;
                case 2:
                    fx = 0; fy = 0; fz = -dz;
                    nx = dx; ny = dy; nz = -dz;
                    return;
                // Need to check up to here for forced!
                case 3:
                    fx = 0; fy = -dy; fz = -dz;
                    nx = dx; ny = -dy; nz = -dz;
                    return;
                case 4:
                    fx = -dx; fy = 0; fz = -dz;
                    nx = -dx; ny = dy; nz = -dz;
                    return;
                case 5:
                    fx = -dx; fy = -dy; fz = 0;
                    nx = -dx; ny = -dy; nz = dz;
                    return;
                // Extras
                case 6:
                    fx = -dx; fy = 0; fz = 0;
                    nx = -dx; ny = 0; nz = dz;
                    return;
                case 7:
                    fx = -dx; fy = 0; fz = 0;
                    nx = -dx; ny = dy; nz = 0;
                    return;
                case 8:
                    fx = 0; fy = -dy; fz = 0;
                    nx = 0; ny = -dy; nz = dz;
                    return;
                case 9:
                    fx = 0; fy = -dy; fz = 0;
                    nx = dx; ny = -dy; nz = 0;
                    return;
                case 10:
                    fx = 0; fy = 0; fz = -dz;
                    nx = 0; ny = dy; nz = -dz;
                    return;
                case 11:
                    fx = 0; fy = 0; fz = -dz;
                    nx = dx; ny = 0; nz = -dz;
                    return;
            }
    }
}