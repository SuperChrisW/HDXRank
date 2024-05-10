solution "Newsim"
configurations { "Release", "Debug" }

newoption {
   trigger     = "tclap",
   value       = "path",
   description = "location of TCLAP (ex: ../tclap-1.2.0)",
}

newoption {
   trigger     = "ticpp",
   value       = "path",
   description = "location of TinyXML++ (ex: ../ticpp-read-only)",
}

if not _OPTIONS["tclap"] then
   _OPTIONS["tclap"] = "../tclap-1.2.0"
end


if not _OPTIONS["ticpp"] then
   _OPTIONS["ticpp"] = "../ticpp-read-only"
end

tclap = _OPTIONS["tclap"]
ticpp = _OPTIONS["ticpp"]
tclapinclude=path.join(tclap,"include")
ticpplibdir=path.join(ticpp,"lib")

configuration { "Release" }
	defines { "NDEBUG" }
	flags   { "OptimizeSpeed" }

configuration { "Debug" }
	defines { "DEBUG" }
	flags { "Symbols" }


project "core"
  language "C++"
  kind     "StaticLib"
  location "lib"
  targetdir "lib"
  files  { "*.cpp" } 
  excludes { "main.cpp" }
  includedirs { tclapinclude, ticpp }
  libdirs { ticpplibdir } 
  links { "m" }
  configuration { "Debug" }
	links { "ticppd" }
  configuration { "Release" }
	links { "ticpp" }

project "main"
  language "C++"
  kind     "ConsoleApp"
  files  { "main.cpp" } 
  includedirs { ticpp }
  links { "m", "core" }
  libdirs { ticpplibdir }
  configuration { "Debug" }
	links { "ticppd" }
  configuration { "Release" }
	links { "ticpp" }

project "measureRMSF"
  language "C++"
  kind     "ConsoleApp"
  location "util"
  targetdir "util"
  files { "util/measureRMSF.cpp" }
  includedirs { "util", "./", tclapinclude }
  links { "core", "m" }

project "measureDistanceTraveled"
  language "C++"
  kind     "ConsoleApp"
  location "util"
  targetdir "util"
  files { "util/measureDistanceTraveled.cpp" }
  includedirs { "util", "./", tclapinclude }
  links { "core", "m" }

project "assignAtomTypes"
  language "C++"
  kind     "ConsoleApp"
  location "util"
  targetdir "util"
  files { "util/assignAtomTypes.cpp", "util/AtomCategories.cpp" }
  includedirs { "util", "./", tclapinclude }
  links { "core", "m" }

project "matchState1ToState2"
  language "C++"
  kind     "ConsoleApp"
  location "util"
  targetdir "util"
  files { "util/matchState1ToState2.cpp", "util/Swap.cpp", "util/AtomLookup.cpp" }
  includedirs { "util", "./", tclapinclude }
  links { "core", "m" }

project "projections2d"
  language "C++"
  kind     "ConsoleApp"
  location "util"
  targetdir "util"
  files { "util/projections2d.cpp" }
  includedirs { "util", "./", tclapinclude }
  links { "core", "m" }
