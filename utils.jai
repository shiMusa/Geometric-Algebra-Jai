CHECKS :: true;

println :: inline (msg: string, args: ..Any) {
    print(msg, ..args);
    print("\n");
}

println :: inline (arg: Any) {
    print("%\n", arg);
}

NewResizableArray :: (count: s64, $T: Type, $initialized := true, allocator: Allocator = .{}, alignment: s32 = -1) -> [..]T {
    arr := NewArray(count, T, initialized, allocator, alignment); // ! alloc
    res := resizable(arr);
    res.count = 0;
    return res;
}

powi :: (base: int, exponent: int) -> int {
    res := 1;
    for 0..exponent-1
        res *= base;
    return res;
}

hamming_weight :: (num: $T) -> int {
    count := 0;
    while num {
        num &= num-1;
        count += 1;
    }
    return count;
}


factorial :: (n: int) -> int {
    res := 1;
    for 2..n res *= it;
    return res;
}

// ? Scientific Computing, Vol I: Linear and nonlinear equations
// ? Texts in computational science and engineering 18
// ? Springer
// ? 2.3.1.10 Binomial Coefficients, page 55
binomial :: (from: int, choose: int) -> int {
    #if CHECKS {
        assert(choose <= from, "Out of range: from = %, choose = %", from, choose);
    }

    if choose == 0 || from == choose {
        return 1;
    }

    return binomial(from - 1, choose) + binomial(from - 1, choose - 1);
}


get_const_value :: (type: Type, name: string, $T: Type) -> (T, bool) {
    s := cast(*Type_Info_Struct)type;
    itype := type_info(T);
    for s.specified_parameters {
        // println("get_const_value: it = %, it.type.type = %, itype.type = %", it, it.type.type, itype.type);
        if it.type.type != itype.type || it.name != name then continue;

        data_segment := s.constant_storage_buffer;
        offset := it.offset_into_constant_storage;
        assert(offset >= 0);

        pointer := cast(*T)(data_segment + offset);
        return <<pointer, true;
    }
    return 0, false;
}

set_const_value :: (type: *Type_Info_Struct, name: string, value: $T) {
    itype := type_info(T);
    for type.specified_parameters {
        if it.type.type != itype.type || it.name != name then continue;

        data_segment := type.constant_storage_buffer;
        offset := it.offset_into_constant_storage;
        assert(offset >= 0);

        pointer := cast(*T)(data_segment + offset);
        <<pointer = value;
    }
}