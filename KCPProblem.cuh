#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iterator>
#include <string>


// Point Structure, used in each Dataset P,Q
typedef struct point //
{
	unsigned long long id;		// Point ID
	double x;			        // Dimension x
    double y;                   // Dimension y
    //unsigned int feature;	// Feature 


} point;

typedef struct pair{

    int startP;
    int startQ;
    int splitP;
    int splitQ;

} pair;


typedef struct closestpairs{

    point p;
    point q;
    double dist;

} closestpairs;

typedef struct ovpairs{

    int startP;
    int startQ;
    int splitP;
    int splitQ;
    int priority;

} ovpairs;


bool endsWith(const std::string& str, const std::string& suffix)
{
    return str.size() >= suffix.size() &&
        str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

typedef std::vector<point> point_vector_t;
typedef std::vector<pair> pair_vector_t;


/** \brief Overloaded operator used for reading point from a text file
 *
 * \param i istream& input stream
 * \param p Point& point to read
 * \return istream& input stream
 *
 */
 std::istream& operator >>(std::istream& i, point& p)
 {
     i >> p.id;
     i >> p.x;
     i >> p.y;
 
     return i;
 }


/** \brief Class definition for kcp problem
 */
class KCPProblem{

public:
    KCPProblem(const std::string& fileDatasetP, const std::string& fileDatasetQ, size_t nofkcp)
        : DatasetP(new point_vector_t), DatasetQ(new point_vector_t){
        
        //set the filenames and read the data files
        this->fileDatasetP = fileDatasetP;
        this->fileDatasetQ = fileDatasetQ;
        this->nofkcp = nofkcp;
        
        this->LoadDataFiles();
    }

    virtual ~KCPProblem(){
    }

    const point_vector_t& GetDatasetP() const{
        return *DatasetP;
    }

    const point_vector_t& GetDatasetQ() const{
        return *DatasetQ;
    }

    size_t GetNofKCP() const{
        return nofkcp;
    }

    virtual size_t GetDatasetPSize() const{
        return DatasetP->size();
    }

    virtual size_t GetDatasetQSize() const{
        return DatasetQ->size();
    }

protected:
    std::string fileDatasetP;
    std::string fileDatasetQ;
    //std::chrono::duration<double> loadingTime;

    /** \brief Template method for loading data files. It can add the points in an internal or external memory vector
     *
     * \param filename const string& the filename to read from
     * \param dataset PointVector& the internal or external memory vector to add points into
     *
     */
    template<class PointVector>
    void LoadFile(const std::string& filename, PointVector& dataset)
    {
        //Handle both binary and text dataset files
        if (endsWith(filename, ".bin"))
        {
            LoadBinaryFile(filename, dataset);
        }
        else if (endsWith(filename, ".txt"))
        {
            LoadTextFile(filename, dataset);
        }
        else{
            LoadCsvFile(filename, dataset);
        }
    }

private:

    size_t nofkcp = 0;
    std::unique_ptr<point_vector_t> DatasetP;
    std::unique_ptr<point_vector_t> DatasetQ;

    void LoadDataFiles(){

        //Record the time for loading the data files
        //auto start = std::chrono::high_resolution_clock::now();

        LoadFile(fileDatasetP, *DatasetP);
        LoadFile(fileDatasetQ, *DatasetQ);

        //auto finish = std::chrono::high_resolution_clock::now();
        //loadingTime = finish - start;
    }

    template<class PointVector>
    void LoadBinaryFile(const std::string& filename, PointVector& dataset)
    {
        //open binary file
        std::fstream fs(filename, std::ios::in | std::ios::binary);
        size_t numPoints = 0;
        //read the number of points at the beginning of the file
        fs.read(reinterpret_cast<char*>(&numPoints), std::streamsize(sizeof(size_t)));
        dataset.reserve(numPoints);

        //read each point and add to vector
        for (size_t i = 0; i < numPoints && !fs.eof(); ++i)
        {
            point p;
            fs.read(reinterpret_cast<char*>(&p), std::streamsize(sizeof(point)));
            dataset.push_back(p);
        }

        fs.close();
    }

    template<class PointVector>
    void LoadTextFile(const std::string& filename, PointVector& dataset)
    {
        std::fstream fs(filename, std::ios::in);
        //std::copy(std::istream_iterator<point>(fs), std::istream_iterator<point>(), std::back_inserter(dataset));
        std::string str;
        int count=0;

        std::getline(fs, str, fs.widen(','));

        while(!fs.eof()){
            
            point newPoint;

            newPoint.id = count;
            newPoint.x = atof(str.c_str());
            std::getline(fs, str, fs.widen('\n'));
            newPoint.y = atof(str.c_str());
            
            dataset.push_back(newPoint);

            count++;
            std::getline(fs, str, fs.widen(','));

        }
        
        
        fs.close();
    }

    template<class PointVector>
    void LoadCsvFile(const std::string& filename, PointVector& dataset)
    {
        std::fstream fs(filename, std::ios::in);
        //std::copy(std::istream_iterator<point>(fs), std::istream_iterator<point>(), std::back_inserter(dataset));
        std::string str;
        int count=0;

        std::getline(fs, str, fs.widen('\t'));

        while(!fs.eof()){
            
            point newPoint;

            newPoint.id = atof(str.c_str());
            std::getline(fs, str, fs.widen('\t'));
            newPoint.x = atof(str.c_str());
            std::getline(fs, str, fs.widen('\n'));
            newPoint.y = atof(str.c_str());
            
            dataset.push_back(newPoint);

            count++;
            std::getline(fs, str, fs.widen('\t'));

        }
        
        
        fs.close();
    }

};
