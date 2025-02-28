#include <AMReX_EB2.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <chrono>
#include <thread>

#include "../EBGeometry/EBGeometry.hpp"

using namespace amrex;

// Template for the AMReX SDF interface using EBGeometry
template <class T, class Meta, class BV, size_t K>
class STLGeometryIF
{
public:
    // Constructor from STL file
    STLGeometryIF(const std::string& filename, bool use_bvh, bool reverse_normals,
                  const Vector<Real>& translate, Real scale, const Vector<Real>& center)
    {
        Print() << "Loading STL file: " << filename << "\n";
        
        try {
            // Load the STL file into a DCEL mesh
            Print() << "Reading STL into DCEL mesh...\n";
            auto mesh = EBGeometry::Parser::readIntoDCEL<T, Meta>(filename);
            Print() << "DCEL mesh created, face count: " << mesh->getFaces().size() << "\n";
            
            // Apply transformations as needed
            if (reverse_normals) {
                Print() << "Reversing STL normals\n";
                mesh->flip();
            }
            
            // Create the SDF representation based on user preference
            if (use_bvh) {
                Print() << "Creating BVH acceleration structure...\n";
                m_sdf = EBGeometry::Parser::readIntoLinearBVH<T, Meta, BV, K>(filename);
                Print() << "BVH acceleration structure created\n";
            } else {
                Print() << "Using standard mesh representation\n";
                m_sdf = EBGeometry::Parser::readIntoMesh<T, Meta>(filename);
            }
            
            // Store transformation parameters for later use in operator()
            m_translate = translate;
            m_scale = scale;
            m_center = center;
            
            // Test the SDF with a simple query to make sure it works
            using Vec3 = EBGeometry::Vec3T<T>;
            Vec3 testPoint(0.0, 0.0, 0.0);
            Print() << "Testing SDF with point (0,0,0), value = " << m_sdf->value(testPoint) << "\n";
            
            Print() << "STL geometry loaded successfully\n";
        }
        catch (const std::exception& ex) {
            Abort("Error loading STL geometry: " + std::string(ex.what()));
        }
    }
    
    // AMReX's implicit function interface
    Real operator()(AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
        using Vec3 = EBGeometry::Vec3T<T>;
        
        // Debug counter - print every 100000 calls to detect hanging
        static int query_count = 0;
        if (++query_count % 100000 == 0) {
            Print() << "SDF query #" << query_count << " at point (" 
                   << x << ", " << y << ", " << z << ")\n";
        }
        
        // Add a safety timeout in case of infinite recursion
        static auto last_output = std::chrono::steady_clock::now();
        auto now = std::chrono::steady_clock::now();
        if (std::chrono::duration_cast<std::chrono::seconds>(now - last_output).count() > 5) {
            Print() << "SDF evaluation taking too long, point = ("
                   << x << ", " << y << ", " << z << ")\n";
            last_output = now;
        }
        
        // Apply transformations before querying the SDF
        Real transformed_x = (x - m_center[0]) * m_scale + m_center[0] + m_translate[0];
        Real transformed_y = (y - m_center[1]) * m_scale + m_center[1] + m_translate[1];
        Real transformed_z = (z - m_center[2]) * m_scale + m_center[2] + m_translate[2];
        
        try {
            return m_sdf->value(Vec3(transformed_x, transformed_y, transformed_z));
        } catch (const std::exception& e) {
            Print() << "Exception in SDF evaluation: " << e.what() << "\n";
            return 1.0; // Return a safe value
        }
    }
    
    // Alternative interface using RealArray
    inline Real operator()(const RealArray& p) const noexcept
    {
        return this->operator()(AMREX_D_DECL(p[0], p[1], p[2]));
    }

private:
    // SDF representation
    std::shared_ptr<EBGeometry::ImplicitFunction<T>> m_sdf;
    
    // Transformation parameters
    Vector<Real> m_translate;
    Real m_scale;
    Vector<Real> m_center;
};

// Main initialization function
void initialize_EB2(const Geometry& geom, const int required_level, const int max_coarsening_level)
{
    BL_PROFILE("initialize_EB2");
    
    // Parse EB settings from input file
    ParmParse pp("eb2");
    std::string geom_type;
    pp.get("geom_type", geom_type);
    
    if (geom_type == "user_defined") {
        // Get STL file parameters
        bool use_bvh = true;
        std::string stl_file = "Shapes.stl";
        int num_coarsen_opt = 0;
        Vector<Real> stl_translation(AMREX_SPACEDIM, 0.0);
        Real stl_scale = 1.0;
        Vector<Real> stl_center(AMREX_SPACEDIM, 0.0);
        bool stl_reverse_normal = false;
        
        pp.query("use_bvh", use_bvh);
        pp.query("stl_file", stl_file);
        pp.query("num_coarsen_opt", num_coarsen_opt);
        pp.queryarr("stl_translation", stl_translation);
        pp.query("stl_scale", stl_scale);
        pp.queryarr("stl_center", stl_center);
        pp.query("stl_reverse_normal", stl_reverse_normal);
        
        Print() << "\n===== STL Geometry Parameters =====\n";
        Print() << "File: " << stl_file << "\n";
        Print() << "Using BVH: " << (use_bvh ? "yes" : "no") << "\n";
        Print() << "Scale: " << stl_scale << "\n";
        Print() << "Center: (";
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            Print() << stl_center[i] << (i < AMREX_SPACEDIM-1 ? ", " : "");
        }
        Print() << ")\n";
        
        Print() << "Translation: (";
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            Print() << stl_translation[i] << (i < AMREX_SPACEDIM-1 ? ", " : "");
        }
        Print() << ")\n";
        
        Print() << "Reverse Normals: " << (stl_reverse_normal ? "yes" : "no") << "\n";
        Print() << "====================================\n\n";
        
        // Set up the STL geometry
        constexpr int K = 4;  // Tree degree
        using T = float;      // Precision for EBGeometry
        using Meta = EBGeometry::DCEL::DefaultMetaData;
        using BV = EBGeometry::BoundingVolumes::AABBT<T>;
        
        // Create the geometry object
        Print() << "Creating STL geometry object...\n";
        STLGeometryIF<T, Meta, BV, K> stl_geometry(stl_file, use_bvh, stl_reverse_normal, 
                                                  stl_translation, stl_scale, stl_center);
        
        // Create the shop and build the EB2 geometry
        Print() << "Creating EB2 shop...\n";
        auto gshop = EB2::makeShop(stl_geometry);
        
        Print() << "Starting EB2::Build...\n";
        // Added extra debug parameters to show progress and reduce work
        EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 4, 
                  true, // verbose
                  true, // max_grid_size
                  num_coarsen_opt);
        
        Print() << "EB Geometry initialization complete\n";
    }
    else {
        // Fall back to regular AMReX initialization for other geometry types
        EB2::Build(geom, max_coarsening_level, max_coarsening_level, 4);
    }
}
