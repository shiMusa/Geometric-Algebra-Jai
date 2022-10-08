const std = @import("std");
const print = std.debug.print;
const ga = @import("geometric_algebra.zig");

const GA = ga.Algebra(3, 0, 1, f32);
const e0 = GA.e(0);
const e1 = GA.e(1);
const e2 = GA.e(2);
const e4 = GA.e(4);

fn Foo(comptime N: u8, comptime B: [N]u8) type {
    return struct {
        pub const N = N;
        pub const bases = B;
        vals: [N]i32,
    };
}

pub fn main() !void {
    const N = 5;
    const b = [N]u8{ 0, 1, 2, 3, 4 };
    var buffer = [N]i32{ 0, -42, 313, 9000, 2022 };
    const F = Foo(N, b);
    var f = F{ .vals = buffer };
    print("f = {}\n", .{f});
    print("F: {}, {any}; sizeOf = {}\n", .{ F.N, F.bases, @sizeOf(F) });

    print("\nGeometric Algebra\n\n", .{});
    print("algebra = {}, B = {any}\n", .{ GA, GA.B });
    print("BN0 = {}\n", .{e0});

    const c = e1.geo(e2);
    print("c = {}\n", .{c});
    const d = e1.geo(e1);
    print("d = {}\n", .{d});
    const e = e4.geo(e4);
    print("e = {}\n", .{e});
}
