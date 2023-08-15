#include "amr.h"

#include "config.h"

#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <PltFileManager.H>

#include <array>

struct AMRState {
    Config        config;
    amrex::AMReX* amr;

    AMRState(Config c) : config(c) {
        // dummy params
        std::string exe_name = "amr2vdb";
        std::array<char*, 1> argv_src = { exe_name.data() };
        int argc = argv_src.size();
        char** argv = argv_src.data();

        amr = amrex::Initialize(argc, argv);
    }

    bool init(std::filesystem::path path) {
        if (!std::filesystem::exists(path)) {
            std::cerr << "Path does not exist: " << path << std::endl;
            return false;
        }

        pele::physics::pltfilemanager::PltFileManager plt_data(path);
        auto plt_vars = plt_data.getVariableList();
        int  nvars    = plt_vars.size();

        std::set<std::string> extant_names(plt_vars.begin(), plt_vars.end());

        for (auto const& vname : config.variables) {
            if (!extant_names.contains(vname)) {
                std::cerr << "Variable " << vname
                          << " does not exist in pltfile.\n";
                return false;
            }
        }

        if (config.max_level < 0) { config.max_level = plt_data.getNlev() - 1; }

        if (config.max_level >= plt_data.getNlev()) {
            std::cerr << "Asking for refinement level " << config.max_level
                      << " but pltfile only has " << plt_data.getNlev()
                      << std::endl;
            return false;
        }

        auto geom = plt_data.getGeom(config.max_level);

        auto domain = geom.Domain();
        auto larray = domain.smallEnd().toArray();
        auto harray = domain.bigEnd().toArray();

        std::cout << "Level:    " << config.max_level << std::endl;
        std::cout << "Domain L: " << larray[0] << " " << larray[1] << " "
                  << larray[2] << std::endl;
        std::cout << "Domain H: " << harray[0] << " " << harray[1] << " "
                  << harray[2] << std::endl;

        auto grid = amrex::BoxArray(domain);

        auto const dmap = amrex::DistributionMapping(grid);

        auto output = amrex::MultiFab(grid, dmap, nvars, 0);

        return true;
    }

    ~AMRState() {
        amrex::Finalize(amr);
    }
};

std::shared_ptr<AMRState> load_file(std::filesystem::path path,
                                    Config const&         c) {
    auto ret = std::make_shared<AMRState>(c);

    if (ret->init(path)) {
        return ret;
    }

    return {};
}

/*

std::shared_ptr<AMRState> load_file(std::filesystem::path) {
amrex::Initialize(argc,argv);
{
  ParmParse pp;
  std::string file;
  pp.get("file",file);
  DataServices::SetBatchMode();
  Amrvis::FileType fileType(Amrvis::NEWPLT);
  DataServices dataServices(file, fileType);
  if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
  }
  AmrData& amrData = dataServices.AmrDataRef();

  int finestLevel = amrData.FinestLevel();
  pp.query("finestLevel",finestLevel);

  const std::vector<std::string>& pieces = Tokenize(file,std::string("/"));

  int slicedir, sliceloc;
  pp.get("slicedir",slicedir);
  pp.get("sliceloc",sliceloc);
  std::string varname;
  pp.get("varname",varname);
  AMREX_ALWAYS_ASSERT(amrData.CanDerive(varname));

  Box domain = amrData.ProbDomain()[finestLevel];
  domain.setSmall(slicedir,sliceloc);
  domain.setBig(slicedir,sliceloc);
  FArrayBox data(domain,1);
  amrData.FillVar(&data,data.box(),finestLevel,varname,0);
  Print() << "min,max: " << data.min() << ", " << data.max() << std::endl;

  std::string outtype("image"); pp.query("outtype",outtype);
  if (outtype=="image" || outtype=="gray")
  {
      Real data_min = data.min();
      Real data_max = data.max();
      pp.query("min",data_min);
      pp.query("max",data_max);

      BaseFab<int> image;
      const int nVals = 256;
      pixelizeData(data,slicedir,sliceloc,image,data_min,data_max,nVals);

      const int width = image.box().length(0);
      const int height= image.box().length(1);

      if (outtype=="image")
      {
          std::string palfile;
          pp.get("palette",palfile);
          int r[256], g[256], b[256], t[256];
          LOAD_PALETTE_STR(palfile,r,g,b,t);

          std::string outfile = pieces[pieces.size()-1] + PPM;
          pp.query("outfile",outfile);
          STORE_PPM_STR(outfile, width, height, image.dataPtr(), r, g, b);
      }
      else
      {
          std::string outfile = pieces[pieces.size()-1] + PGM;
          pp.query("outfile",outfile);
          STORE_PGM_STR(outfile, width, height, image.dataPtr());
      }
  }
}
amrex::Finalize();
}
*/
