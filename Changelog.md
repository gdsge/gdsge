# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.2] - Unreleased

### Added

- Editor plugins

### Changed

- Add instructions for compiling a gmod file in README
- Add instructions for setting up compilers on a macOS
- Fix include dependence bugs
- Fix bugs in using var_tensor

## [0.1.1] - 2022-09-30

### Changed

- Remove dependencies on v2struct in convert_to_interp_eval_array.m. This should improve performance

## [0.1.0] - 2022-09-09

### Added

- Changelog
- cinclude <> 
- Blocks with snippets added to mex_xxx.cpp: start_loop, finish_loop
- Blocks with snippets added to iter_xxx.cpp: pre_mex_call, pre_jac_code, post_jac_code

### Changed

- Ending of cxx; block is changed to endcxx;

### Removed

- Figure in Mendoza example that produces broken titles
