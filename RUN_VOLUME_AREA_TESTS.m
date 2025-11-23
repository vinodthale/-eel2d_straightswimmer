%% ================================================================
%  RUN ALL VOLUME AND AREA TESTS
%
%  This script runs comprehensive volume and area verification tests
%  for the eel2d geometry generation.
%
%  Tests include:
%  1. Main eel geometry generation with area/volume calculations
%  2. Detailed volume and area verification test
%% ================================================================

clear all; clc; close all;

fprintf('\n');
fprintf('================================================================\n');
fprintf('  EEL2D VOLUME AND AREA TEST SUITE\n');
fprintf('================================================================\n\n');

%% ================================================================
%  TEST 1: Main Eel Geometry with Area/Volume Calculations
%% ================================================================
fprintf('Running Test 1: Main eel geometry generation...\n');
fprintf('----------------------------------------------------------------\n\n');

try
    run('eel2d_straightswimmer.m');
    fprintf('\n✓ Test 1 PASSED: Eel geometry generated successfully\n\n');
catch ME
    fprintf('\n✗ Test 1 FAILED: %s\n\n', ME.message);
    rethrow(ME);
end

%% ================================================================
%  TEST 2: Detailed Volume and Area Verification
%% ================================================================
fprintf('================================================================\n');
fprintf('Running Test 2: Detailed volume and area verification...\n');
fprintf('----------------------------------------------------------------\n\n');

try
    run('test_volume_area.m');
    fprintf('\n✓ Test 2 PASSED: Volume and area verification completed\n\n');
catch ME
    fprintf('\n✗ Test 2 FAILED: %s\n\n', ME.message);
    rethrow(ME);
end

%% ================================================================
%  SUMMARY
%% ================================================================
fprintf('================================================================\n');
fprintf('  TEST SUITE SUMMARY\n');
fprintf('================================================================\n\n');

fprintf('All tests completed successfully!\n\n');

fprintf('Generated files:\n');
fprintf('  ├── output/vertices/eel2d.vertex\n');
fprintf('  ├── output/images/eel2d_straightswimmer.png\n');
fprintf('  ├── output/images/volume_area_verification.png\n');
fprintf('  ├── output/data/eel_geometry_properties.txt\n');
fprintf('  └── output/data/volume_area_data.txt\n\n');

fprintf('Check the output directory for detailed results.\n\n');

fprintf('================================================================\n');
fprintf('  TEST SUITE COMPLETE\n');
fprintf('================================================================\n\n');
