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


// * mergesort for counting the number of neighbor-swaps necessary to sort an array
// * see https://www.geeksforgeeks.org/number-swaps-sort-adjacent-swapping-allowed/

count_swaps :: (arr: []$B) -> int {
    tmp := NewArray(arr.count, B, false);
    defer array_free(tmp);

    return merge_sort(arr, tmp, 0, arr.count-1);
}

merge_sort :: (arr: []$B, tmp: []B, left: int, right: int) -> int {
    inv_count := 0;
    if right > left {
        mid := (right + left)/2;

        inv_count  = merge_sort(arr, tmp, left, mid);
        inv_count += merge_sort(arr, tmp, mid+1, right);

        inv_count += merge(arr, tmp, left, mid+1, right);
    }
    return inv_count;
}

merge :: (arr: []$B, tmp: []B, left: int, mid: int, right: int) -> int {
    inv_count := 0;

    i := left;
    j := mid;
    k := left;
    while i <= mid-1 && j <= right {
        if arr[i] <= arr[j] {
            tmp[k] = arr[i];
            k += 1;
            i += 1;
        } else {
            tmp[k] = arr[j];
            k += 1;
            j += 1;
            inv_count += mid - i;
        }
    }

    while i <= mid-1 {
        tmp[k] = arr[i];
        k += 1;
        i += 1;
    }

    while j <= right {
        tmp[k] = arr[j];
        k += 1;
        j += 1;
    }

    for i: left..right {
        arr[i] = tmp[i];
    }

    return inv_count;
}


basis_to_array :: (b: $B) -> []int {
    res : [..]int;
    x := b;
    i := 0;
    while x {
        if x & 1 {
            array_add(*res, i);
        }
        i += 1;
        x = x >> 1;
    }
    return res;
}

count_swaps :: (a: $B, b: B) -> int {
    DEBUG :: false;

    #if DEBUG println("count swaps of % and %", a, b);

    arr_a := basis_to_array(a);
    defer array_free(arr_a);
    if arr_a.count == 0 then return 0;

    arr_b := basis_to_array(b);
    defer array_free(arr_b);
    if arr_b.count == 0 then return 0;

    combined := NewArray(arr_a.count + arr_b.count, int, false);
    defer array_free(combined);

    for arr_a combined[              it_index] = it;
    for arr_b combined[arr_a.count + it_index] = it;
    #if DEBUG println("combined = %", combined);

    n_swaps := count_swaps(combined);
    #if DEBUG println("n_swaps = %", n_swaps);

    return n_swaps;
}