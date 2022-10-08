const std = @import("std");
const print = std.debug.print;

pub fn Algebra(comptime POS: u8, comptime NEG: u8, comptime ZER: u8, comptime T: type) type {
    return struct {
        pub const N = POS + NEG + ZER;
        pub const POS = POS;
        pub const NEG = NEG;
        pub const ZER = ZER;
        pub const T = T;
        pub const B = switch (N) {
            0...8 => u8,
            9...16 => u16,
            17...32 => u32,
            else => u64,
        };

        pub const one: T = 1;
        pub const zero: T = 0;

        fn basis_to_T(a: B, b: B) bool {
            return (a & b) >> (POS + NEG) > 0 or a ^ b == 0;
        }

        fn geo_type(a: B, b: B) type {
            if (basis_to_T(a, b)) return T;
            return BasisNumber(a | b);
        }

        fn e_basis(comptime dim: B) B {
            if (dim == 0) return 0;
            return 1 << dim - 1;
        }
        pub fn e(comptime dim: B) BasisNumber(e_basis(dim)) {
            return BasisNumber(e_basis(dim)){ .value = one };
        }

        pub fn BasisNumber(comptime basis: B) type {
            return struct {
                pub const basis = basis;
                value: T,

                pub fn geo(self: @This(), b: anytype) geo_type(basis, @TypeOf(b).basis) {
                    // const _B = @TypeOf(b);
                    const RES = geo_type(basis, @TypeOf(b).basis);
                    if (RES == T) {
                        if (basis >> (POS + NEG) > 0) return zero; // basis squares to 0
                        return self.value * b.value;
                    }
                    return RES{ .value = self.value * b.value };
                }
            };
        }
    };
}
