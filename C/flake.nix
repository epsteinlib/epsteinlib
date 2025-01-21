# SPDX-FileCopyrightText: 2024 Jan Schmitz <schmitz@num.uni-sb.de>
#
# SPDX-License-Identifier: AGPL-3.0-only
{
  outputs = {
    self,
    nixpkgs,
    flake-parts,
    pre-commit-hooks-nix,
    devshell,
    ...
  } @ inputs: let
    supportedSystems = [
      "x86_64-linux"
      "aarch64-linux"
      "x86_64-darwin"
      "aarch64-darwin"
    ];
    nixpkgsConfig = {
      allowUnsupportedSystem = true;
    };
    flake_parts_out = flake-parts.lib.mkFlake {inherit inputs;} {
      imports = [
        inputs.pre-commit-hooks-nix.flakeModule
        inputs.treefmt-nix.flakeModule
        inputs.devshell.flakeModule
        ./nix/packages.nix
        ./nix/devshells.nix
      ];
      systems = supportedSystems;
      perSystem = {
        self',
        config,
        inputs',
        lib,
        pkgs,
        system,
        ...
      }: rec {
        _module.args.pkgs = import nixpkgs {
          inherit system;
          config = nixpkgsConfig;
        };

        devShells = let
          inherit (self'.devShells) _epstein;
        in rec {
          epstein = pkgs.mkShell {
            inherit (_epstein) shellHook;

            inputsFrom = [
              self'.packages.default
            ];
            packages = [
              self'.packages.default
              self'.packages.default.optional-dependencies.dev
            ];
          };
          epstein_devshell = pkgs.mkShell {
            shellHook = epstein.shellHook;

            inputsFrom = [
              epstein
            ];
            packages = [
              config.pre-commit.settings.enabledPackages
            ];
          };
          default = epstein_devshell;
        };

        treefmt = {
          projectRootFile = "flake.nix";
          programs = {
            # JSON
            biome.enable = true;
            # C
            clang-format.enable = true;
            # Nix
            alejandra.enable = true;
            # Python
            black.enable = true;
            # TOML
            taplo.enable = true;
            # Yaml
            yamlfmt.enable = true;
          };
        };
      };
    };
  in
    flake_parts_out;

  inputs = {
    nixpkgs-stable.url = "github:nixos/nixpkgs/nixos-24.05";
    nixpkgs-unstable.url = "github:nixos/nixpkgs/nixos-unstable";
    nixpkgs.follows = "nixpkgs-stable";

    flake-parts.url = "github:hercules-ci/flake-parts";

    pre-commit-hooks-nix.url = "github:cachix/git-hooks.nix";

    treefmt-nix = {
      url = "github:numtide/treefmt-nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };

    devshell = {
      url = "github:numtide/devshell";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };
}
