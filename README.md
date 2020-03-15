### Usage

    root -l -b -q main_proc.C+'("input_file.list","output_file.root")'
    
where `"input_file.list"` - filelist of mcpico data, `"output_file.root"` - resulted file.
Resolution calculation:

    root -l -q -b res2.C
    
In `res2.C` one should specify `"output_file.root"` to calculate resolution correction factor from.
