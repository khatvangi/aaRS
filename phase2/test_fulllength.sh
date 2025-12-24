#!/bin/bash
cd inputs/af3_jsons_fulllength
/storage/kiran-stuff/alphafold3/af3_wrapper_optimized.sh --json_path=fulllength_shallow_pro.json 2>&1 | tee ../../logs/test_fulllength.log
