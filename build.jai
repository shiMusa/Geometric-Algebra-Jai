#import "Basic";
#import "Compiler";

#run build();

build :: () {
    w := compiler_create_workspace("Target Program");
    if !w {
        print("Workspace creation failed.\n");
        return;
    }
    
    options := get_build_options(w);

    args := options.compile_time_command_line;
    print("\nargs: %\n", args);
    filename := args[1];

    options.output_executable_name = filename;
    
    // * activate inlining
    options.enable_bytecode_inliner = true;

    // options.optimization_level = .RELEASE;
    
    set_build_options(options, w);

    compiler_begin_intercept(w);
    add_build_file(filename, w);
    message_loop();
    compiler_end_intercept(w);

    set_build_options_dc(.{do_output=false});
}

message_loop :: () {
    while true {
        message := compiler_wait_for_message();

        if message.kind == {
            case .COMPLETE;
                break;
        }
    }
}