{
  description = "A development flake for GreenX. Usage `nix develop .#devShell.x86_64-linux`";
  nixConfig.bash-prompt = ''\[\033[1;31m\][\[\033[0m\]\[\033[1;37m\]dev .\[\033[0m\]\[\033[1;34m\] $(basename \$$PWD)\[\033[0m\]\[\033[1;31m\]]\[\033[0m\] '';

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs = { self, nixpkgs }:
    let
      pkgs = nixpkgs.legacyPackages.x86_64-linux;

      zofu = pkgs.stdenv.mkDerivation {
        pname = "zofu";
        version = "1.1.1";

        src = pkgs.fetchFromGitHub {
          owner = "acroucher";
          repo = "zofu";
          rev = "96cf65829610c7812f212b42268081af356c8f52";
          sha256 = "1i9ylxy52x7gg95dvrdb72pb45l2rsnrlnsi504cvpyngfnnxd0n";
        };

        buildInputs = [ pkgs.cmake pkgs.gfortran ];

        configurePhase = ''
          mkdir build
          cd build
          cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$out ../
        '';

        buildPhase = ''
          make -j4
        '';

        installPhase = ''
          make install
        '';
        
        meta = {
          homepage = "https://github.com/acroucher/zofu";
          description = "A Fortran unit test framework.";
        };
      };

      pygreenx = pkgs.python3Packages.buildPythonPackage {
        pname = "pygreenx";
        version = "1.0";

        src = ./python;

        propagatedBuildInputs = with pkgs.python3Packages; [
          pytest
          setuptools
          numpy
        ];

        meta = {
          homepage = "https://github.com/nomad-coe/greenX";
          description = "Python utils for the GreenX library.";
        };
      };

      custom_python_env = pkgs.python3.withPackages (ps: [
        pygreenx
      ]);
    in
      {
        devShell.x86_64-linux = pkgs.mkShell {
          buildInputs = with pkgs; [
            # Python packages
            custom_python_env
            # Additional libraries and utils
            zofu
            python310Packages.pytest
            doxygen
            lapack-reference
            cmake
            valgrind
            gmp
          ];
        };
      };
}
