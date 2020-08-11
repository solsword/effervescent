"use strict";

// ---------
// Constants
// ---------

// Global canvas context
var CTX = undefined;

// Minimum zoom-in
var MIN_SCALE = 8;

// Maximum zoom-out
var MAX_SCALE = 64;

// Maximum automatic zoom-out
var AUTO_MAX = 32;

// Comfortable scale level
var COMFORTABLE_SCALE = MIN_SCALE * 2;

// How much bigger than the interest bounding box should the scale be
var IDEAL_SCALE_MULTIPLIER = 1.6;

// Speed at which to change scales (percentage of scale difference per second)
var ZOOM_IN_SPEED = 1.3
var ZOOM_OUT_SPEED = 2.5

// Speed at which to pan the origin (percentage of distance-to-ideal-origin per
// second)
var PAN_SPEED = 1.2;

// The default color of the maze
var DEFAULT_MAZE_COLOR = "#cccccc";

// Colors to use for successive distance-from-origin values
var DISTANCE_COLORS = [
    "#d088ff",

    "#cc88ff",
    "#c888ff",
    "#c488ff",
    "#c088ff",
    "#bc88ff",
    "#b888ff",
    "#b488ff",
    "#b088ff",
    "#ac88ff",
    "#a888ff",
    "#a488ff",
    "#a088ff",
    "#9c88ff",
    "#9888ff",
    "#9488ff",
    "#9088ff",
    "#8c88ff",

    "#8888ff",

    "#888cff",
    "#8890ff",
    "#8894ff",
    "#8898ff",
    "#889cff",
    "#88a0ff",
    "#88a4ff",
    "#88a8ff",
    "#88acff",
    "#88b0ff",
    "#88b4ff",
    "#88b8ff",
    "#88bcff",
    "#88c0ff",
    "#88c4ff",
    "#88c8ff",
    "#88ccff",

    "#88d0ff",

    "#88ccff",
    "#88c8ff",
    "#88c4ff",
    "#88c0ff",
    "#88bcff",
    "#88b8ff",
    "#88b4ff",
    "#88b0ff",
    "#88acff",
    "#88a8ff",
    "#88a4ff",
    "#88a0ff",
    "#889cff",
    "#8898ff",
    "#8894ff",
    "#8890ff",
    "#888cff",

    "#8888ff",

    "#8c88ff",
    "#9088ff",
    "#9488ff",
    "#9888ff",
    "#9c88ff",
    "#a088ff",
    "#a488ff",
    "#a888ff",
    "#ac88ff",
    "#b088ff",
    "#b488ff",
    "#b888ff",
    "#bc88ff",
    "#c088ff",
    "#c488ff",
    "#c888ff",
    "#cc88ff",
];

// The color of the destination
var DEST_COLOR = "#ffffff";

// The color of the destination
var POS_COLOR = "#3399ff";

// The palette
var PALETTE = [
  "#ff4444", // red
  "#ffff22", // yellow
  "#4466ff", // blue
  "#ff44cc", // pink
  "#44ccff", // light aqua
  "#ffaa22", // orange
  "#66ff66", // green
  "#f8ffaa", // cream
  "#bbeeff", // sky blue
  "#aaff44", // lime green
];

// The seed for the maze:
var MAZE_SEED = 1947912873;

// Size of each pattern
var PATTERN_SIZE = 2;

// Length of each loop
var LOOP_LENGTH = (PATTERN_SIZE * 4) - 4;

// Object placeholder for things we're in the process of generating
var WORKING_ON_IT = {};

// The fractally increasing layer caches. Values in the process of being
// generated are represented by WORKING_ON_IT, while never-requested
// values are undefined.
var LAYER_CACHES = {};

// The per-seed caches for distance-from-origin information
var DIST_CACHES = {};

// Queue for layers waiting to be generated. Each entry should be a set
// of fractal coordinates. Coordinates in the queue which cannot be
// generated because the superstructure they belong to hasn't been
// created yet will be removed and discarded.
var GEN_QUEUE = [];

// Queue for positions waiting to have their distance computed. Each
// entry is a triple containing an age counter, a seed, and a set of
// absolute coordinates.
var DIST_QUEUE = [];

// Number of distances to process during a single distance computation
// step.
var DIST_BATCH = 12;

// How many steps before we throw old unsuccessful entries out of the
// DIST_QUEUE.
var DIST_AGE_LIMIT = 100;

// Number of layers to generate per gen step.
var GEN_SPEED = 12;

// Delay (ms) between generation ticks
var GEN_DELAY = 2;

// Delay (ms) between test attempts
var TEST_DELAY = 50;

// Delay (ms) between automatic steps
var STEP_DELAY = 20;

// Minimum and maximum callback delays (in ms):
var MIN_DELAY = 1;
var MAX_DELAY = 2000;

// How long to wait between auto destination checks
var AUTO_DEST_DELAY = 30;

// How long to wait (in AUTO_DEST_DELAY cycles) before automatically setting a
// new destination.
var AUTO_DEST_WAIT = 100;

// Counter to keep track of the number of cycles until we should automatically
// set a new random destination.
var AUTO_DEST_COUNTER = 0;

// Orientations
var NORTH = 1;
var EAST = 2;
var SOUTH = 4;
var WEST = 8;

// Which frame we're on:
var FRAME = 0;

// When the frame counter resets:
var MAX_FC = 1000;

// Have we reached the destination?
var AT_DESTINATION = true;

// Current mode
var MODE = undefined;

// -------------------------
// Updaters & Event Handlers
// -------------------------

function set_mode(mode) {
    MODE = mode;
}

function set_scale(context, scale_factor) {
    // Scale is in world-units-per-canvas-width
    if (scale_factor < MIN_SCALE) {
        scale_factor = MIN_SCALE;
    }
    let alt_max = MAX_SCALE * context.cwidth / context.cheight;
    let limit = Math.min(MAX_SCALE, alt_max);
    if (scale_factor > limit) {
        scale_factor = limit;
    }
    context.scale = scale_factor;
}

function set_origin(context, origin) {
    context.origin = origin;
}

function set_position(context, grid_coords) {
    // Sets the current position
    context.position = grid_coords;
}

function set_destination(context, grid_coords) {
    // Sets the current destination.
    context.destination = grid_coords;
}

function update_canvas_size(canvas, context) {
    // Updates the canvas size. Called on resize after a timeout.
    var bounds = canvas.getBoundingClientRect();
    var car = bounds.width / bounds.height;
    canvas.width = 800 * car;
    canvas.height = 800;
    context.cwidth = canvas.width;
    context.cheight = canvas.height;
    context.middle = [context.cwidth / 2, context.cheight / 2];
    context.bounds = bounds;
}

// Scrolling constants
var PIXELS_PER_LINE = 18;
var LINES_PER_PAGE = 40;

function handle_scroll(ctx, ev) {
    let unit = ev.deltaMode;
    let dy = ev.deltaY;

    // Normalize units to pixels:
    if (unit == 1) {
        dy *= PIXELS_PER_LINE;
    } else if (unit == 2) {
        dy *= PIXELS_PER_LINE * LINES_PER_PAGE;
    }

    set_scale(ctx, ctx.scale * Math.max(0.1, (100 + dy) / 100));
}

function event_pos(ctx, ev) {
    // Returns viewport position of event.
    if (ev.touches) {
        ev = ev.touches[0];
    }
    return pgc__vc(ctx, [ev.clientX, ev.clientY]);
}

function on_canvas(vc) {
    return (
        0 <= vc[0] && vc[0] <= 1
        && 0 <= vc[1] && vc[1] <= 1
    );
}

function handle_tap(ctx, ev) {
    let vc = event_pos(ctx, ev);
    if (on_canvas(vc)) {
        let gc = wc__gc(cc__wc(ctx, vc__cc(ctx, vc)));
        set_destination(ctx, gc);
    }
}

// --------------------
// Conversion functions
// --------------------

// Page <-> viewport coordinates
function pgc__vc(ctx, pc) {
    return [
        (pc[0] - ctx.bounds.left) / ctx.bounds.width,
        (pc[1] - ctx.bounds.top) / ctx.bounds.height
    ];
}

function vc__pgc(ctx, vc) {
    return [
        ctx.bounds.left + ctx.bounds.width * vc[0],
        ctx.bounds.top + ctx.bounds.height * vc[1],
    ];
}

// Viewport <-> canvas coordinates
function vc__cc(ctx, vc) {
    return [
        vc[0] * ctx.cwidth,
        vc[1] * ctx.cheight
    ];
}

function cc__vc(ctx, cc) {
    return [
        cc[0] / ctx.cwidth,
        cc[1] / ctx.cheight
    ];
}

// Canvas <-> world coordinates
function cc__wc(ctx, cc) {
    return [
        ((cc[0] - ctx.cwidth/2)/ctx.cwidth) * ctx.scale + ctx.origin[0],
        -((cc[1] - ctx.cheight/2)/ctx.cwidth) * ctx.scale + ctx.origin[1]
            // scale ignores canvas height
    ];
}

function wc__cc(ctx, wc) {
    return [
        ((wc[0] - ctx.origin[0]) / ctx.scale) * ctx.cwidth + ctx.cwidth/2,
        -((wc[1] - ctx.origin[1]) / ctx.scale) * ctx.cwidth + ctx.cheight/2
    ];
}

function canvas_unit(ctx) {
    // Returns the length of one world-coordinate unit in canvas coordinates.
    return (ctx.cwidth / ctx.scale);
}

// World <-> grid coordinates
function wc__gc(wc) {
    return [
        Math.floor(wc[0] + 0.5),
        Math.floor(wc[1] + 0.5)
    ];
}

function gc__wc(gc) {
    return [
        gc[0],
        gc[1]
    ];
}

// Page coordinates all the way to grid coordinates:
function pgc__gc(ctx, pgc) {
    return wc__gc(
        cc__wc(
            ctx,
            vc__cc(
                ctx,
                pgc__vc(ctx, pgc)
            )
        )
    );
}

// Gets extrema of canvas in the grid. Returns an object with keys 'NW', 'NE',
// 'SW', and 'SE' for each of the four corners.
function grid_extrema(ctx) {
    return {
        'NW': pgc__gc(ctx, [ ctx.bounds.left, ctx.bounds.top ]),
        'NE': pgc__gc(ctx, [ ctx.bounds.right, ctx.bounds.top ]),
        'SW': pgc__gc(ctx, [ ctx.bounds.left, ctx.bounds.bottom ]),
        'SE': pgc__gc(ctx, [ ctx.bounds.right, ctx.bounds.bottom ]),
    };
}

// Edge + socket <-> edge ID
// Undefined socket value indicates unconstrained socket
function ec__eid(ec) {
    if (ec[1] == undefined) {
        return "" + ec[0] + ":A"
    } else {
        return "" + ec[0] + ":" + ec[1];
    }
}

function eid__ec(eid) {
    let ori = parseInt(eid[0]);
    let slot = parseInt(eid.slice(2));
    if (isNaN(slot)) {
        return [ ori, undefined ];
    } else {
        return [ ori, slot ];
    }
}

// (Possibly-)Nonspecific edge coordinate to (possible-singular) list of EIDs.
function ec__eids(ec) {
    let result = [];
    if (ec[1] != undefined) {
        result.push(ec__eid(ec));
    } else {
        for (let i = 0; i < PATTERN_SIZE; ++i) {
            result.push(ec__eid([ec[0], i]));
        }
    }
    return result;
}

// Edge + socket <-> pattern index
// Pattern indices are:
//   0 1
//   3 2
function ec__pidx(ec) {
    if (ec[1] == undefined) {
        console.warn(
            "Attempted to translate an underspecified edge coordiante into a "
          + "pattern coordinate."
        );
        return undefined;
    } else if (ec[0] == NORTH) {
        return ec[1];
    } else if (ec[0] == EAST) {
        return 1 + ec[1];
    } else if (ec[0] == SOUTH) {
        return 3 - ec[1];
    } else { // ec[0] == WEST, we hope
        return (4 - ec[1]) % 4;
    }
}

// horizontal determines whether we return horizontal or vertical edge
// coordinates.
function pidx__ec(pidx, horizontal) {
    if (pidx == 0) {
        if (horizontal) {
            return [ NORTH, 0 ];
        } else {
            return [ WEST, 0 ];
        }
    } else if (pidx == 1) {
        if (horizontal) {
            return [ NORTH, 1 ];
        } else {
            return [ EAST, 0 ];
        }
    } else if (pidx == 2) {
        if (horizontal) {
            return [ SOUTH, 1 ];
        } else {
            return [ EAST, 1 ];
        }
    } else { // pidx == 3, we hope
        if (horizontal) {
            return [ SOUTH, 0 ];
        } else {
            return [ WEST, 1 ];
        }
    }
}

// Pattern absolute x/y <-> pattern index
// Pattern indices are:
//   0 1
//   3 2
// Absolute coordinates in a pattern are:
//   (0, 1)   (1, 1)
//   (0, 0)   (1, 0)
function pidx__pc(pidx) {
    let y = pidx < 2 ? 1 : 0;
    let x = pidx % 3 == 0 ? 0 : 1;
    return [x, y];
}

function pc__pidx(pc) {
    if (pc[0] == 0 && pc[1] == 0) {
        return 3;
    } else if (pc[0] == 1 && pc[1] == 0) {
        return 2;
    } else if (pc[0] == 1 && pc[1] == 1) {
        return 1;
    } else if (pc[0] == 0 && pc[1] == 1) {
        return 0;
    }
}


function ori__vec(ori) {
    // Converts an orientation into an [x, y] absolute coordinate direction
    // vector.
    if (ori == NORTH) {
        return [0, 1];
    } else if (ori == EAST) {
        return [1, 0];
    } else if (ori == SOUTH) {
        return [0, -1];
    } else if (ori == WEST) {
        return [-1, 0];
    } else {
        console.error("Bad orientation: " + ori);
    }
}

function edge_ori(edge) {
    // Returns the direction in which an edge extends from the associated
    // edge absolute coordinates, which will be either SOUTH (for EAST
    // and WEST edges) or EAST (for NORTH and SOUTH edges).
    if (edge == EAST || edge == WEST) {
        return SOUTH;
    } else {
        return EAST;
    }
}

// ------------
// Drawing Code
// ------------

function draw_frame(now) {
    // Draws a single frame & loops itself

    // Measure time
    let ms_time = window.performance.now();
    if (CTX.now == undefined) {
        CTX.now = ms_time;
        window.requestAnimationFrame(draw_frame);
        return; // skip this frame to get timing for the next one
    }
    CTX.elapsed = ms_time - CTX.now;
    CTX.now = ms_time;

    // Count frames
    FRAME += 1;
    FRAME %= MAX_FC;

    // Clear the canvas:
    CTX.clearRect(0, 0, CTX.cwidth, CTX.cheight);

    adjust_viewport(CTX);
    draw_maze(CTX, MAZE_SEED);
    draw_destination(CTX);
    draw_position(CTX);

    // Requeue ourselves
    if (!FAILED) {
        window.requestAnimationFrame(draw_frame);
    } else {
        console.error("Draw loop aborted due to test failure.");
    }
}

function interest_bb(ctx) {
    // Computes the bounding box of the interesting region (the region
    // containing all points of each trail) in world coordinates.
    let result = {
        "left": ctx.destination[0],
        "right": ctx.destination[0],
        "top": ctx.destination[1],
        "bottom": ctx.destination[1]
    }

    let pos = ctx.position;
    if (pos[0] < result.left) { result.left = pos[0]; }
    if (pos[0] > result.right) { result.right = pos[0]; }
    if (pos[1] < result.top) { result.top = pos[1]; }
    if (pos[1] > result.bottom) { result.bottom = pos[1]; }

    return result;
}

function adjust_viewport(ctx) {
    // Adjusts the scaling factor and origin according to points of
    // interest
    let ibb = interest_bb(ctx);

    // Compute ideal scale
    let ar = (ctx.cwidth / ctx.cheight);
    let ideal_scale = Math.max(
        COMFORTABLE_SCALE,
        ibb.right - ibb.left,
        (ibb.bottom - ibb.top) * ar,
    ) * IDEAL_SCALE_MULTIPLIER;

    let use_scale = ideal_scale;
    let alt_max = AUTO_MAX * ctx.cwidth / ctx.cheight;
    let scale_limit = Math.min(AUTO_MAX, alt_max);
    if (use_scale > scale_limit) { use_scale = scale_limit; }

    // Compute scale difference
    let scale_diff = use_scale - ctx.scale;

    let zs;
    if (scale_diff > 0) { // zooming out
        zs = ZOOM_OUT_SPEED;
    } else {
        zs = ZOOM_IN_SPEED;
    }

    // Adjust the context scale
    if (MODE == "attract") {
        ctx.scale = Math.max(
            MIN_SCALE,
            ctx.scale + zs * scale_diff * (ctx.elapsed / 1000)
        );
    } else {
        // We don't change the scale
        use_scale = ctx.scale;
    }

    // Compute ideal center
    let ideal_center = [
        (ibb.left + ibb.right) / 2,
        (ibb.top + ibb.bottom) / 2
    ];

    let use_center = ideal_center;
    if (use_scale < ideal_scale) {
        // Adjust center towards the current position to prioritize its
        // visibility.
        let ratio = use_scale / ideal_scale;
        let aim_vector = [
            ideal_center[0] - ctx.position[0],
            ideal_center[1] - ctx.position[1]
        ];
        use_center = [
            ctx.position[0] + aim_vector[0] * ratio,
            ctx.position[1] + aim_vector[1] * ratio
        ];
    }

    let center_diff = [
        use_center[0] - ctx.origin[0],
        use_center[1] - ctx.origin[1]
    ];

    // Pan towards the ideal center
    ctx.origin = [
        ctx.origin[0]
      + PAN_SPEED * center_diff[0] * (ctx.elapsed / 1000),
        ctx.origin[1]
      + PAN_SPEED * center_diff[1] * (ctx.elapsed / 1000)
    ];
}

function draw_maze(ctx, seed) {
    // Draws the visible portion of the labyrinth.

    // Set stroke color:
    ctx.strokeStyle = DEFAULT_MAZE_COLOR;
    ctx.lineWidth = 2;

    // Radius of each grid cell
    let cell_size = canvas_unit(ctx);

    // Iterate over visible (and a few invisible) cells at the base layer:
    let extrema = grid_extrema(ctx);

    for (let x = extrema['SW'][0] - 1; x <= extrema['SE'][0] + 1; ++x) {
        for (let y = extrema['SW'][1] - 1; y <= extrema['NW'][1] + 1; ++y) {
            // Absolute position
            let ac = [x, y];

            // Canvas coordinates for this grid cell:
            let cc = wc__cc(ctx, gc__wc(ac));

            // Draw edge links for each cell
            let fc = ac__fc(seed, ac);
            let connections = full_adjacency_at(fc);
            let distance = lookup_distance_from_origin(seed, ac);
            if (distance == undefined) {
                ctx.strokeStyle = DEFAULT_MAZE_COLOR;
            } else {
                ctx.strokeStyle = DISTANCE_COLORS[
                    distance % DISTANCE_COLORS.length
                ];
            }

            if (connections == undefined) {
                // not available yet; has been requested

                // Just draw a circle
                ctx.beginPath();
                ctx.arc(cc[0], cc[1], cell_size*0.2, 0, 2*Math.PI);
                ctx.stroke();

            } else { // From-links for each direction
                // Assemble a list of neighbor coordinates to connect to
                let neighbors = [];
                if (is_connected(connections, NORTH)) {
                    neighbors.push([ cc[0], cc[1] - cell_size/2 ]);
                }

                if (is_connected(connections, EAST)) {
                    neighbors.push([ cc[0] + cell_size/2, cc[1] ]);
                }

                if (is_connected(connections, SOUTH)) {
                    neighbors.push([ cc[0], cc[1] + cell_size/2 ]);
                }

                if (is_connected(connections, WEST)) {
                    neighbors.push([ cc[0] - cell_size/2, cc[1] ]);
                }

                // Draw simple lines:
                for (let nb of neighbors) {
                    ctx.beginPath();
                    ctx.moveTo(cc[0], cc[1]);
                    ctx.lineTo(nb[0], nb[1]);
                    ctx.stroke();
                }
            }
        }
    }
}

function draw_destination(ctx) {
    // Draws the current destination for the given context.
    ctx.lineWidth = 2 + canvas_unit(ctx) * 0.05;
    ctx.strokeStyle = DEST_COLOR;
    ctx.fillStyle = DEST_COLOR;
    let α = 1;
    if (MODE == 'attract') {
        α = 1 - (AUTO_DEST_COUNTER / AUTO_DEST_WAIT);
    }
    ctx.beginPath();
    let cc = wc__cc(ctx, ctx.destination);
    ctx.arc(cc[0], cc[1], canvas_unit(ctx)*0.2, 0, 2*Math.PI);
    ctx.globalAlpha = α;
    ctx.stroke();
    ctx.globalAlpha = 0.25 * α;
    ctx.fill();
    ctx.globalAlpha = 1;
}

function draw_position(ctx) {
    // Draws the current position for the given context.
    ctx.lineWidth = 2 + canvas_unit(ctx) * 0.05;
    ctx.strokeStyle = POS_COLOR;
    ctx.fillStyle = POS_COLOR;
    ctx.beginPath();
    let cc = wc__cc(ctx, ctx.position);
    ctx.arc(cc[0], cc[1], canvas_unit(ctx)*0.15, 0, 2*Math.PI);
    ctx.globalAlpha = 1;
    ctx.stroke();
    ctx.globalAlpha = 0.25;
    ctx.fill();
    ctx.globalAlpha = 1;
}


// ---------
// Misc Code
// ---------

function blend_color(c1, c2, r) {
    // Blends 1-r of the first color with r of the second color. Does silly RGB
    // interpolation.
    c1 = c1.slice(1);
    c2 = c2.slice(1);

    let r1 = parseInt(c1.slice(0, 2), 16);
    let g1 = parseInt(c1.slice(2, 4), 16);
    let b1 = parseInt(c1.slice(4, 6), 16);
    let a1 = 255;
    if (c1.length > 6) { let a1 = parseInt(c1.slice(6, 8)); }

    let r2 = parseInt(c2.slice(0, 2), 16);
    let g2 = parseInt(c2.slice(2, 4), 16);
    let b2 = parseInt(c2.slice(4, 6), 16);
    let a2 = 255;
    if (c2.length > 6) { let a2 = parseInt(c2.slice(6, 8)); }

    let new_r = Math.floor(r1 * (1 - r) + r2 * r);
    let new_g = Math.floor(g1 * (1 - r) + g2 * r);
    let new_b = Math.floor(b1 * (1 - r) + b2 * r);
    let new_a = Math.floor(a1 * (1 - r) + a2 * r);

    let hr = new_r.toString(16);
    if (hr.length == 1) { hr = "0" + hr; }
    let hg = new_g.toString(16);
    if (hg.length == 1) { hg = "0" + hg; }
    let hb = new_b.toString(16);
    if (hb.length == 1) { hb = "0" + hb; }
    let ha = new_a.toString(16);
    if (ha.length == 1) { ha = "0" + ha; }
    return "#" + hr + hg + hb + ha;
}

function lfsr(x) {
    // Implements a max-cycle-length 32-bit linear-feedback-shift-register.
    // See: https://en.wikipedia.org/wiki/Linear-feedback_shift_register
    // Note that this is NOT reversible!
    var lsb = x & 1;
    var r = x >>> 1;
    if (lsb) {
        r ^= 0x80200003; // 32, 22, 2, 1
    }
    return r;
}

function randint(up_to, seed) {
    // Picks a random integer strictly less than the given value using the given
    // seed.
    return posmod(lfsr(seed), up_to);
}

var FLIP_RESOLUTION = 10000000;

function flip_biased(p, seed) {
    // Flips a biased coin, returning true with probability p (determined
    // by the given seed).
    return randint(FLIP_RESOLUTION, seed) < FLIP_RESOLUTION * p
}

function choose_randomly(possibilities, seed) {
    // Picks randomly from a list using the given seed.
    let idx = posmod(seed, possibilities.length);
    return possibilities[idx];
}

function posmod(n, base) {
    // Mod operator that always returns positive results.
    return ((n % base) + base) % base;
}

// -------------------
// Fractal Coordinates
// -------------------

function origin_for(seed) {
  // Returns the origin coordinates for the given seed. These will always be
  // within 2*PATTERN_SIZE units of the true origin.
  let r = lfsr(seed + 17371947103);
  let x = r % (2*PATTERN_SIZE);
  r = lfsr(r);
  let y = r % (2*PATTERN_SIZE);
  return [ x, y ];
}

// Computes the absolute coordinates of the southwest corner of the
// origin fractal grid unit that's at the given height.
// Fractal grid units are nested as alternately the southwest and
// northeast portions of their parent units, like so:
//
//                                *--*--*
//                                |  |  |
//  height = 0                    *--*--*
//  origin = (0, 0)               |<>|  |
//                                *--*--*
//                         
//                         
//                          *-----*-----*
//                          |     |     |
//                          |     |  *  |
//                          |     |     |
//  height = 1              *-----*-----*
//  origin = (-2, -2)       |     |     |
//                          |     |     |
//                          |     |     |
//                          *-----*-----*
//                         
//                         
//                          *-----------*-----------*
//                          |           |           |
//                          |           |           |
//                          |           |           |
//                          |           |           |
//                          |           |           |
//                          |           |           |
//                          |           |           |
//  height = 2              *-----------*-----------*
//  origin = (-2, -2)       |           |           |
//                          |           |           |
//                          |           |           |
//                          |     *     |           |
//                          |           |           |
//                          |           |           |
//                          |           |           |
//                          *-----------*-----------*
//
function fractal_origin(height) {
    // With a 0, 1, 1, 2, 2, 3, 3, ... pattern, adjacent heights share
    // origins...
    let grouped = Math.ceil(height/2);
    // The offset is the sum of 2, 8, 32, 128, ...
    let offset = 0;
    for (let i = 1; i <= grouped; ++i) {
        offset += Math.pow(PATTERN_SIZE, i*2 - 1);
    }
    return [ -offset, -offset ];
}

// Fractal <-> absolute coordinates
//
// Absolute coordinates are an [x, y] pair that denotes a grid cell on a
// standard right-handed grid where +x goes East and +y goes North.
//
// For subspecific fractal coordinates, returns the absolute coordinates
// of the southwest corner of the cell. Conversion takes time
// proportional to the log of the absolute coordinate's distance from the
// origin.
//
// Fractal coordinates consist of a height value and a list of pattern
// indices indicating a trace downwards from that height. A height of 0
// indicates the unit is an n×n region, 1 an n^2×n^2 region, and so on.
// Each entry in a trace is an index between 0 and 3 that indicates which
// sub-cell of the current cell the location is within. The trace may be
// shorter than the height, in which case the fractal coordinates denote
// a layer above the base grid cells.

function fc__ac(fc) {
    let seed = fc_seed(fc);
    let height = fc_height(fc); // 3
    let trace = fc_trace(fc); // [2, 0]

    let cw = Math.pow(PATTERN_SIZE, height); // width of a single cell
    // cw = 8
    let result = [0, 0]; // bottom-left default

    // Trace down through each layer:
    for (let i = 0; i < trace.length; ++i) {
        let pc = pidx__pc(trace[i]);
        let x = pc[0];
        let y = pc[1];
        result[0] += x * cw;
        result[1] += y * cw;
        cw /= PATTERN_SIZE;
    }

    let origin = origin_for(seed);
    let fo = fractal_origin(height);

    return [
        result[0] + origin[0] + fo[0],
        result[1] + origin[1] + fo[1]
    ];
}


function ac__fc(seed, ac) {
    // Account for seed-specific origin
    let origin = origin_for(seed);
    let off_ac = [ ac[0] - origin[0], ac[1] - origin[1] ];

    let sw = [0, 0];
    let ne = [1, 1];
    let height = 0;
    let scale = 1;
    // Keep going until we're in-bounds
    while (
        off_ac[0] < sw[0]
        || off_ac[0] > ne[0]
        || off_ac[1] < sw[1]
        || off_ac[1] > ne[1]
    ) {
        height += 1;
        scale *= PATTERN_SIZE;
        if (height % 2 == 1) { // An odd step, so we expand to the bottom-left
            sw[0] -= scale;
            sw[1] -= scale;
        } else {
            ne[0] += scale;
            ne[1] += scale;
        }
    }

    // Compute coordinates within the layer at the given height
    let in_ac = [
        off_ac[0] - sw[0],
        off_ac[1] - sw[1]
    ];

    let trace = [];
    for (let h = height; h >= 0; h -= 1) {
        let x = in_ac[0] < scale ? 0 : 1;
        let y = in_ac[1] < scale ? 0 : 1;
        let idx = pc__pidx([x, y]);
        trace.push(idx);
        in_ac[0] -= x * scale;
        in_ac[1] -= y * scale;
        scale /= PATTERN_SIZE;
    }

    return build_fc(seed, height, trace);
}

function fc__edge_ac(fc, edge) {
    // Works like fc__ac, but computes the coordinates of the start of a
    // specific edge of the given layer. The edge is specified using one
    // of the constants NORTH, EAST, SOUTH, or WEST.
    //
    // The coordinates returned will be at the top (north end) of WEST
    // and EAST edges, and at the left (west) end of NORTH and SOUTH
    // edges, and they will identify the north/western cell of the two
    // cells that are adjacent across the edge at that point.

    // Compute height of the fc and height that the edge is at based on
    // trace length. For a fully-specific fc, the edge_height will be 0.
    let height = fc_height(fc);
    let trace = fc_trace(fc);
    let edge_height = fc_layers_below(fc);

    // edge_length for height = 0 will be 1...
    let edge_length = Math.pow(PATTERN_SIZE, edge_height);

    // Southwest corner as identified by fc__ac
    let result = fc__ac(fc);

    // Move from southwest corner to appropriate corner as dictated by
    // the edge paramter (no movement when edge_height is 0 because
    // edge_length = 1). Also flip across the edge for NORTH and WEST
    // edges so that results are the same as for matching SOUTH and EAST
    // edges.
    let across = edge_length - 1;
    if (edge == NORTH) {
        result[1] += across; // across to NW corner
        result[1] += 1; // north one to match SOUTH edge above
    } else if (edge == EAST) {
        result[0] += across; // across to SE corner
        result[1] += across; // across to NE corner
        // WEST edges will move to join us
    } else if (edge == SOUTH) {
        // stay where we are; NORTH edges will join us here
    } else if (edge == WEST) {
        result[1] += across; // across to NW corner
        result[0] -= 1; // west one to match EAST edge on our left
    } else {
        console.error("Invalid edge value in fc__edge_ac: " + edge);
    }

    return result;
}

function fc_seed(fc) {
    // Returns the seed of the given fractal coordinates
    return fc[0];
}

function fc_height(fc) {
    // Returns the height of the given fractal coordinates
    return fc[1];
}

function fc_trace(fc) {
    // Returns the trace of the given fractal coordinates.
    return fc[2];
}

function fc_layers_below(fc) {
    // Computes the number of unspecified layers for the given fractal
    // coordinates, which is just the heigh minus the length of the
    // trace, plus 1. A coordinate with height 1 and a single trace entry
    // has 0 layers below it, one with height 2 and a single trace entry
    // has 1 layer below it, etc.
    return fc_height(fc) + 1 - fc_trace(fc).length;
}

function build_fc(seed, height, trace) {
    // Builds a set of fractal coordinates out of a seed, a height, and a
    // trace.
    return [seed, height, trace];
}

function clone_fc(fc) {
    // Returns a clone of the given fractal coordinates.
    return build_fc(
        fc_seed(fc),
        fc_height(fc),
        fc_trace(fc).slice()
    );
}

function build_child_fc(parent_fc, pidx) {
    // Builds the fractal coordinates for the child of the given
    // coordinates with the given pattern index.
    return build_fc(
        fc_seed(parent_fc),
        fc_height(parent_fc),
        fc_trace(parent_fc).concat([pidx])
    );
}

function extend_fc(fc) {
    // Extends the given fractal coordinates so that their height is
    // increased by one while still denoting the same cell. Returns a new
    // set of coordinates without modifying the originals.
    let seed = fc_seed(fc);
    let height = fc_height(fc);
    let trace = fc_trace(fc);
    let new_idx;
    if (height % 2 == 0) {
        // even steps encapsulate as northeast corner of new layer that
        // extends farther southwest.
        new_idx = 1;
    } else {
        // odd steps encapsulate as southwest corner of new layer that extends
        // farther northeast.
        new_idx = 3;
    }
    return [
        seed,
        height + 1,
        [ new_idx ].concat(trace)
    ];
}

function parent_of(fc) {
    // Returns the fractal coordinates of the parent of the given cell.
    let seed = fc_seed(fc);
    let height = fc_height(fc);
    let trace = fc_trace(fc);
    if (trace.length <= 1) {
        return parent_of(extend_fc(fc));
    }
    return [
        seed,
        height,
        trace.slice(0, trace.length - 1)
    ];
}

function normalize_fc(fc) {
    // Retracts the given fractal coordinates as much as possible,
    // getting rid of unnecessary height and central indices. The result
    // still refers to the same cell. Returns a new set of coordinates
    // without modifying the originals. The result will have a trace
    // length of at least 1, even where it could technically have a trace
    // length of 0, unless it is given a fractal coordinate with trace
    // length 0, which it will return as-is.
    let seed = fc_seed(fc);
    let height = fc_height(fc);
    let trace = fc_trace(fc).slice();

    if (height == 0 || trace.length <= 1) {
        // Cannot be shortened
        return fc;
    }

    // If we chop off the top & then call extend_fc, what index will
    // extend_fc generate?
    let natural_inner_idx;
    if (height % 2 == 1) {
        natural_inner_idx = 1;
    } else {
        natural_inner_idx = 3;
    }

    if (trace[0] == natural_inner_idx) {
        // can be shortened; recurse
        return normalize_fc(build_fc(seed, height - 1, trace.slice(1)));
    } else {
        // cannot be further shortened
        return fc;
    }
}

function local_seed(fr_coords) {
    // Determines the local layer seed for the given fractal coordinate
    // location.
    let seed = fc_seed(fr_coords);
    let height = fc_height(fr_coords);
    let trace = fc_trace(fr_coords);
    for (let i = 0; i < height; ++i) {
        seed = lfsr(seed);
    }
    for (let pidx of trace) {
        seed = lfsr(seed + (seed+1)*pidx);
    }

    return seed;
}

function edge_seed(fr_coords, edge) {
    // Determines the seed for the given edge of the fractally specified
    // layer. Will return the same seed for both fractal + edge
    // coordinates that reference each edge. For example,
    //
    // The EAST edge of [ 192892, 1, [0, 1] ]
    // ...and the WEST edge of [ 192892, 1, [1, 0] ]
    //
    // are the same edge. Similarly:
    //
    // The NORTH edge of [ 38928, 2, [2] ]
    // ...and the SOUTH edge of [ 38928, 2, [1] ]
    //
    // are the same edge, even though they aren't at height 0.
    //
    let seed = fc_seed(fr_coords);
    let ec = fc__edge_ac(fr_coords, edge);
    let edge_height = fc_layers_below(fr_coords) + 1;
    let mix = ((seed + (17*ec[0])) ^ ec[1]) + 3*edge_height;
    let churn = 1 + mix % 4;
    for (let i = 0; i < churn; ++i) {
        mix = lfsr(mix + edge_height);
    }
    return mix;
}

function natural_lower_coords(seed, height) {
    // Returns the fractal coordinates for the natural layer just below
    // the given height (or just the southwest grid cell for height=0).
    // The returned fractal coordinates have a height equal to the given
    // height.
    let result = build_fc(seed, 0, [3]);
    while (fc_height(result) < height) {
        result = extend_fc(result);

    }
    return build_fc(seed, height, fc_trace(result).slice(0, 1));
}

// Generate an array of destinations near the origin
var CLOSE_DESTINATIONS = [];
var CLOSE_SIZE = 20;
for (let x = -CLOSE_SIZE; x <= CLOSE_SIZE; ++x) {
    for (let y = -CLOSE_SIZE; y <= CLOSE_SIZE; ++y) {
        CLOSE_DESTINATIONS.push([x, y]);
    }
}
// Fisher-Yates to randomize order in which we visit destinations
var RDEST_SEED = lfsr(39872081);
for (let i = 0; i < CLOSE_DESTINATIONS.length; ++i) {
    let remaining = CLOSE_DESTINATIONS.length - i - 1;
    let choice = posmod(RDEST_SEED, remaining + 1);
    let tmp = CLOSE_DESTINATIONS[i];
    CLOSE_DESTINATIONS[i] = CLOSE_DESTINATIONS[i + choice];
    CLOSE_DESTINATIONS[i + choice] = tmp;
    RDEST_SEED = lfsr(RDEST_SEED + lfsr(i));
}

// Which nearby destination we're going to now.
var WHICH_DEST = 0;

function next_destination() {
    // Returns the next destination from a shuffled list of nearby
    // destinations.
    WHICH_DEST = posmod(WHICH_DEST + 1, CLOSE_DESTINATIONS.length);
    return CLOSE_DESTINATIONS[WHICH_DEST];
}

// ------------------
// Caching and Lookup
// ------------------

function request_central_layer(seed) {
    if (!LAYER_CACHES.hasOwnProperty(seed)) {
        LAYER_CACHES[seed] = [];
    }
    let cache = LAYER_CACHES[seed];
    if (cache[cache.length - 1] == WORKING_ON_IT) {
        // We're already working on the next central layer
        return;
    }
    let height = cache.length;
    cache[height] = WORKING_ON_IT;
    let fr_coords = natural_lower_coords(seed, height);
    GEN_QUEUE.push(fr_coords);
}

function lookup_layer(fr_coords) {
    // Looks up the cached layer at the given fractal coordinates, or returns
    // undefined and adds an entry to the generation queue if that layer or one
    // of its parents is not yet cached.
    let seed = fc_seed(fr_coords);
    let height = fc_height(fr_coords);
    let trace = fc_trace(fr_coords);

    let cache = LAYER_CACHES[seed];
    if (cache == undefined) {
        cache = [];
        LAYER_CACHES[seed] = cache;
    }

    // If the cache isn't tall enough yet, request that it grow taller
    // and return undefined
    if (cache.length < height + 1) {
        request_central_layer(seed);
        return undefined;
    }
    let ancestor = cache[height];
    if (ancestor == WORKING_ON_IT) {
        return undefined;
    }

    let sofar = [];
    for (let pidx of trace) {
        sofar.push(pidx);
        if (ancestor.children[pidx] == WORKING_ON_IT) {
            // we're already working on it
            return undefined;
        } else if (ancestor.children[pidx] == undefined) {
            // add this to our generation queue
            ancestor.children[pidx] = WORKING_ON_IT;
            GEN_QUEUE.push(build_fc(seed, height, sofar));
            return undefined;
        } else {
            // inwards; onwards
            ancestor = ancestor.children[pidx];
        }
    }
    // we've found it!
    return ancestor;
}

function gen_step() {
    // Self-queuing function that processes the generation queue.
    for (let i = 0; i < GEN_SPEED; ++i) {
        gen_next();
    }
    window.setTimeout(gen_step, GEN_DELAY);
}

function gen_next() {
    // Generates the next layer in the generation queue, first generating a
    // single extra level of the layer cache if at least one more has been
    // requested.
    let next = GEN_QUEUE.shift();
    if (next == undefined) {
        return; // nothing to do right now
    }
    let seed = fc_seed(next);
    let cache = LAYER_CACHES[seed];
    if (cache == undefined) {
        cache = [];
        LAYER_CACHES[seed] = cache;
    }
    if (cache.length < 1 || cache[cache.length - 1] == WORKING_ON_IT) {
        let above;
        above = gen_central_layer(seed, cache.length - 1);
        cache[cache.length - 1] = above;
    }
    let height = fc_height(next);
    let trace = fc_trace(next);
    // Pop last entry in trace (points to layer we're being asked to generate)
    // and keep the rest to find our parent:
    let last = trace[trace.length - 1];
    let remaining = trace.slice(0, trace.length - 1);
    let parent_fc = build_fc(seed, height, remaining);
    let parent = lookup_layer(parent_fc);
    if (
        parent != WORKING_ON_IT
     && parent != undefined
     && parent.children != null
     && parent.children[last] == WORKING_ON_IT
    ) {
        parent.children[last] = gen_layer(parent_fc, parent, last);
    }
}

// --------------------
// Distance Computation
// --------------------

function lookup_distance_from_origin(seed, ac) {
    // Looks up the distance from the origin to the given absolute
    // coordinate position in the grid with the given seed. Returns
    // undefined and requests computation if that distance has not yet
    // been computed.
    let cache = DIST_CACHES[seed];
    if (cache == undefined) {
        cache = {};
        DIST_CACHES[seed] = cache;
    }

    let result = cache["" + ac];
    if (result == undefined) {
        // Queue it for processing if it's novel
        queue_if_novel([0, seed, ac], DIST_QUEUE);
        return undefined;
    } else {
        return result;
    }
}

function queue_if_novel(entry, queue) {
    let already = false;
    for (let existing of queue) {
        if (
            entry[2][0] == existing[2][0]
         && entry[2][1] == existing[2][1]
         && entry[1] == existing[1]
        ) {
            // we reset the age of the match
            existing[0] = 0;
            already = true;
            break;
        }
    }
    if (!already) {
        queue.push(entry);
    }
}

function dist_step() {
    // Self-queuing function that processes the distances queue.
    dist_next(DIST_BATCH);
    window.setTimeout(dist_step, GEN_DELAY);
}

function dist_next(limit) {
    // Computes the distance for the next valid position in the current
    // distances queue, skipping positions whose distances cannot be
    // deduced from already-cached neighbors. Quits after processing
    // limit entries, or when it runs out of entries that can be
    // processed.
    if (DIST_QUEUE.length == 0) {
        return; // nothing to do right now
    }

    let leftovers = [];

    let processed = 0;

    for (let idx = 0; idx < DIST_QUEUE.length; ++idx) {
        let entry = DIST_QUEUE[idx];
        let age = entry[0];
        let seed = entry[1];
        let ac = entry[2];

        let revised_entry = [age + 1, seed, ac];

        if (age > DIST_AGE_LIMIT) {
            // Ignore and drop this entry if it's too old
            continue;
        }

        if (processed >= limit) {
            // no more processing allowed this iteration
            leftovers.push(revised_entry);
            continue;
        }

        let cache = DIST_CACHES[seed];
        if (cache == undefined) {
            cache = {};
            DIST_CACHES[seed] = cache;
        }

        // This *is* the origin!
        if (ac[0] == 0 && ac[1] == 0) {
            cache["" + ac] = 0;
            continue;
        }

        // Gather adjacency information
        let adj = full_adjacency_at(ac__fc(seed, ac));
        if (adj == undefined) {
            // Not enough info to process this entry right now
            leftovers.push(revised_entry);
            continue;
        }

        // Update best-known-distance based on cached values from
        // neighbors
        let best_dist = undefined;
        let neighbors = [];
        for (let dir of [NORTH, EAST, SOUTH, WEST]) {
            if (is_connected(adj, dir)) {
                let nv = ori__vec(dir);
                let nb = [ ac[0] + nv[0], ac[1] + nv[1] ];
                neighbors.push(nb);
                let nbdist = cache["" + nb];
                if (
                    nbdist != undefined
                 && (
                        best_dist == undefined
                     || best_dist > nbdist + 1
                    )
                ) {
                    best_dist = nbdist + 1;
                }
            }
        }

        // Update the cache if appropriate
        let old = cache["" + ac];
        if (best_dist != undefined && (old == undefined || best_dist < old)) {
            // Count this as a processed entry
            processed += 1;
            cache["" + ac] = best_dist;
            // Since we're changing a cached value, re-queue our
            // neighbors for updates.
            for (let nb of neighbors) {
                if (cache["" + nb] != undefined) {
                    queue_if_novel([0, seed, nb], leftovers);
                }
            }
        } else if (best_dist == undefined) {
            // re-queue this entry since we didn't find any useful
            // information from neighbors
            leftovers.push(revised_entry);
        }
    }

    // Replace the queue
    DIST_QUEUE = leftovers;
}


// -------------
// Movement Code
// -------------

function advance_step(ctx) {

    let here_fc = ac__fc(MAZE_SEED, ctx.position);
    let dest_fc = ac__fc(MAZE_SEED, ctx.destination);
    let dir = naive_direction_towards(here_fc, dest_fc);
    if (dir != undefined) {
        if (dir == 0) {
            AT_DESTINATION = true;
        } else {
            AT_DESTINATION = false;
            let vec = ori__vec(dir);
            let next_ac = [
                ctx.position[0] + vec[0],
                ctx.position[1] + vec[1]
            ];
            if (next_ac != undefined) {
                set_position(ctx, next_ac);
            }
        }
    }

    // Requeue
    if (!FAILED) {
        window.setTimeout(advance_step, STEP_DELAY, ctx);
    } else {
        console.error("Stopped destination advance due to test failure.");
    }
}

function adjust_tempo(n) {
    // Adjusts the STEP_DELAY mode by the given multiplier.
    STEP_DELAY *= n;
    if (STEP_DELAY <= MIN_DELAY) { STEP_DELAY = MIN_DELAY; }
    if (STEP_DELAY >= MAX_DELAY) { STEP_DELAY = MAX_DELAY; }
}

function set_tempo(n) {
    // Sets the tempo variable for the given mode to the given value
    // (should be between 0 for min-tempo and 1 for max-tempo. The scale
    // is exponential.
    let escale = 7;
    let x = (Math.exp(escale*n) - 1) / (Math.exp(escale) - 1);
    STEP_DELAY = MIN_DELAY + (MAX_DELAY - MIN_DELAY) * x;
}

function set_tempo_from_slider(slider) {
    set_tempo(1 - parseInt(slider.value) / parseInt(slider.max));
}

function scramble_destination(ctx) {
    if (MODE == "attract" && AT_DESTINATION) {
        AUTO_DEST_COUNTER += 1;
        if (AUTO_DEST_COUNTER >= AUTO_DEST_WAIT) {
            set_destination(ctx, next_destination());
            AUTO_DEST_COUNTER = 0;
        }
    } else {
        AUTO_DEST_COUNTER = 0;
    }

    // Requeue
    if (!FAILED) {
        window.setTimeout(scramble_destination, AUTO_DEST_DELAY, ctx);
    } else {
        console.error("Stopped destination scrambling due to test failure.");
    }
}



function common_parent(from_fc, to_fc) {
    // Finds the fractal coordinates of the closest parent layer which
    // contains both from_fc and to_fc. Returns a list containing those
    // coordinates followed by the pattern indices of the from and to
    // coordinates within that layer. Returns undefined if given
    // coordinates with different seeds. Check seeds:
    if (fc_seed(from_fc) != fc_seed(to_fc)) {
        return undefined;
    }

    // Extend the heights of each coordinate to match:
    while (fc_height(from_fc) < fc_height(to_fc)) {
        from_fc = extend_fc(from_fc);
    }
    while (fc_height(to_fc) < fc_height(from_fc)) {
        to_fc = extend_fc(to_fc);
    }

    // Extend each one more time so they definitely coincide somewhere:
    from_fc = extend_fc(from_fc);
    to_fc = extend_fc(to_fc);

    // Common seed:
    let seed = fc_seed(from_fc);

    // Max-height:
    let height = fc_height(from_fc);

    // Loop to find where they first differ and remember where they're
    // last the same:
    let co_fc = build_fc(seed, height, []);
    let fr_pidx, to_pidx;
    for (let i = 0; i < height + 1; ++i) {
        fr_pidx = fc_trace(from_fc)[i];
        to_pidx = fc_trace(to_fc)[i];
        if (fr_pidx != to_pidx) {
            break;
        } // else extend shared trace:
        fc_trace(co_fc).push(fr_pidx);
    }

    // Return the combined coordinates along with the positions (pattern
    // indices) of the start and end coordinates within that layer.
    return [ normalize_fc(co_fc), fr_pidx, to_pidx ];
}

function is_inside(outer_fc, inner_fc) {
    // Returns whether the given inner_fc is inside the given outer_fc.
    if (fc_seed(outer_fc) != fc_seed(inner_fc)) {
        return false; // seeds don't match
    }

    while (fc_height(outer_fc) < fc_height(inner_fc)) {
        outer_fc = extend_fc(outer_fc);
    }

    let oh = fc_height(outer_fc);
    let ih = fc_height(inner_fc);

    let hd = oh - ih;

    let o_trail = fc_trace(outer_fc);
    let i_trail = fc_trace(inner_fc);

    var dscnt;
    for (dscnt = hd; dscnt < o_trail.length; ++dscnt) {
        if (i_trail[dscnt - hd] != o_trail[dscnt]) {
            return false;
        }
    }
    return dscnt - hd < i_trail.length;
}

function immediate_path_towards_neighbor(fc, pidx, direction) {
    // For the given layer, computes the immediate direction of travel
    // within the sub-layer in order to travel from the sub-layer at the
    // given pattern index out of the layer in the given direction.
    // Returns 0 if there is no outgoing connection from the given layer
    // in the given direction, and undefined if the required information
    // is not yet loaded (that information will be requested
    // automatically).
    //
    // If the given pattern index position is on the edge in the
    // appropriate direction and has an outgoing connection in that
    // direction, then that direction is returned. Otherwise, we return a
    // direction to navigate within the layer to get to a position with
    // those properties.

    let layer = lookup_layer(fc);
    if (layer == undefined) { return undefined; }

    let pattern = layer.pattern;

    let adj = full_adjacency_at(fc);
    if (adj == undefined) { return undefined; }

    if (!is_connected(adj, direction)) {
        // there's no connection in that direction
        return 0;
    }

    let child_fc = build_child_fc(fc, pidx);

    let child_adj = full_adjacency_at(child_fc);
    if (child_adj == undefined) { return undefined; }

    let all_candidates = [];
    let exit_locations = [];

    // Determine which sub-layers are on the given edge
    if (direction == NORTH) {
        all_candidates = [0, 1];
    } else if (direction == EAST) {
        all_candidates = [1, 2];
    } else if (direction == SOUTH) {
        all_candidates = [2, 3];
    } else {
        all_candidates = [3, 0];
    }

    // Determine which sub-layers along the given edge connect out
    for (let exit_pidx of all_candidates) {
        let exit_fc = build_child_fc(fc, exit_pidx);
        let exit_adj = full_adjacency_at(exit_fc);
        if (exit_adj == undefined) { return undefined; }
        if (is_connected(exit_adj, direction)) {
            exit_locations.push(exit_pidx);
        }
    }

    if (exit_locations.length == 0) {
        // No connections in the given direction; should be impossible
        // based on previous check, but who knows?
        return 0;
    }

    // Check if we're already on an exit tile?
    for (let exit_pidx of exit_locations) {
        if (pidx == exit_pidx) {
            // We can immediately exit in the given direction
            return direction;
        }
    }

    // At this point we know we're not on an exit.

    // Local function to test for an exit:
    function is_exit(pidx) {
        return exit_locations.includes(pidx);
    }

    // If we're not on an exit, figure out where to go...
    // Note that this code *assumes* total connectivity within each
    // layer!
    if (direction == NORTH) {
        if (pidx == 0) {
            // on an edge, but not an exit
            if (is_connected(child_adj, EAST)) {
                return EAST;
            } else {
                return SOUTH;
            }
        } else if (pidx == 1) {
            // on an edge, but not an exit
            if (is_connected(child_adj, WEST)) {
                return WEST;
            } else {
                return SOUTH;
            }
        } else if (pidx == 2) {
            // not on edge and not an exit
            if (
                is_connected(child_adj, NORTH)
             && (is_exit(1) || connected_in_pattern(pattern, 1, 0))
            ) {
                return NORTH;
            }
            return WEST; // only other option
        } else if (pidx == 3) {
            // not on edge and not an exit
            if (
                is_connected(child_adj, NORTH)
             && (is_exit(0) || connected_in_pattern(pattern, 0, 1))
            ) {
                return NORTH;
            }
            return EAST; // only other option
        } else {
            console.error("Bad pidx", pidx);
        }
    } else if (direction == EAST) {
        if (pidx == 1) {
            // on an edge, but not an exit
            if (is_connected(child_adj, SOUTH)) {
                return SOUTH;
            } else {
                return WEST;
            }
        } else if (pidx == 2) {
            // on an edge, but not an exit
            if (is_connected(child_adj, NORTH)) {
                return NORTH;
            } else {
                return WEST;
            }
        } else if (pidx == 3) {
            // not on edge and not an exit
            if (
                is_connected(child_adj, EAST)
             && (is_exit(2) || connected_in_pattern(pattern, 2, 1))
            ) {
                return EAST;
            }
            return NORTH; // only other option
        } else if (pidx == 0) {
            // not on edge and not an exit
            if (
                is_connected(child_adj, EAST)
             && (is_exit(1) || connected_in_pattern(pattern, 1, 2))
            ) {
                return EAST;
            }
            return SOUTH; // only other option
        } else {
            console.error("Bad pidx", pidx);
        }
    } else if (direction == SOUTH) {
        if (pidx == 2) {
            // on an edge, but not an exit
            if (is_connected(child_adj, WEST)) {
                return WEST;
            } else {
                return NORTH;
            }
        } else if (pidx == 3) {
            // on an edge, but not an exit
            if (is_connected(child_adj, EAST)) {
                return EAST;
            } else {
                return NORTH;
            }
        } else if (pidx == 0) {
            // not on edge and not an exit
            if (
                is_connected(child_adj, SOUTH)
             && (is_exit(3) || connected_in_pattern(pattern, 3, 2))
            ) {
                return SOUTH;
            }
            return EAST; // only other option
        } else if (pidx == 1) {
            // not on edge and not an exit
            if (
                is_connected(child_adj, SOUTH)
             && (is_exit(2) || connected_in_pattern(pattern, 2, 3))
            ) {
                return SOUTH;
            }
            return WEST; // only other option
        } else {
            console.error("Bad pidx", pidx);
        }
    } else if (direction == WEST) {
        if (pidx == 3) {
            // on an edge, but not an exit
            if (is_connected(child_adj, NORTH)) {
                return NORTH;
            } else {
                return EAST;
            }
        } else if (pidx == 0) {
            // on an edge, but not an exit
            if (is_connected(child_adj, SOUTH)) {
                return SOUTH;
            } else {
                return EAST;
            }
        } else if (pidx == 1) {
            // not on edge and not an exit
            if (
                is_connected(child_adj, WEST)
             && (is_exit(0) || connected_in_pattern(pattern, 0, 3))
            ) {
                return WEST;
            }
            return SOUTH; // only other option
        } else if (pidx == 2) {
            // not on edge and not an exit
            if (
                is_connected(child_adj, WEST)
             && (is_exit(3) || connected_in_pattern(pattern, 3, 0))
            ) {
                return WEST;
            }
            return NORTH; // only other option
        } else {
            console.error("Bad pidx", pidx);
        }
    } else {
        console.error("Bad direction:", direction);
    }

    console.error("Fell out of cases.");
}

function naive_direction_towards(from_fc, to_fc) {
    // Computes the direction (NORTH, EAST, SOUTH, or WEST) that's
    // required to travel from the from_fc to the to_fc (both fractal
    // coordinates). Returns undefined if the required information is not
    // yet loaded. Returns 0 if the two coordinates identify the same
    // absolute position, or if one of the two specified layers is
    // contained within the other.
    //
    // This function reads the fractal connectivity information from the
    // layer data and figures out a path which isn't necessarily the
    // shortest path, as it doesn't take possible-but-unlikely shortcuts
    // through higher-level layers into account, nor does it reason
    // intelligently about alternatives when multiple paths exist.

    // Check for matching coordinates:
    let from_ac = fc__ac(from_fc);
    let to_ac = fc__ac(to_fc);
    if (same(from_ac, to_ac)) {
        return 0;
    }

    // Find common ancestor
    let [common_fc, from_pidx, to_pidx] = common_parent(from_fc, to_fc);

    let common_layer = lookup_layer(common_fc);
    if (common_layer == undefined) {
        return undefined;
    }

    if (from_pidx == to_pidx) {
        // parent/child relationship -> no movement: from_pidx is already
        // inside to_pidx, or vice versa.
        return 0;
    } else {
        // separate ancestries
        // compute direction of movement within pattern
        let dir = direction_in_pattern(
            common_layer.pattern,
            from_pidx,
            to_pidx
        );

        // Set up to iterate downwards until we get to a fully-specific
        // fractal coordinate
        let below = fc_layers_below(common_fc);
        let extended_from_fc = from_fc;
        while (fc_height(extended_from_fc) < fc_height(common_fc)) {
            extended_from_fc = extend_fc(extended_from_fc);
        }
        let from_trace = fc_trace(extended_from_fc);

        // now compute local direction of movement by continuously
        // applying immediate_path_towards_neighbor.
        for (
            let tidx = from_trace.length - below + 1;
            tidx < from_trace.length;
            ++tidx
        ) {
            let specific_fc = build_fc(
                fc_seed(extended_from_fc),
                fc_height(extended_from_fc),
                from_trace.slice(0, tidx) // tidx itself excluded
            );
            dir = immediate_path_towards_neighbor(
                specific_fc,
                from_trace[tidx],
                dir
            );
        }

        return dir;
    }
}

function direction_towards(from_fc, to_fc) {
    // Computes the direction (NORTH, EAST, SOUTH, or WEST) that's
    // required to travel from the from_fc to the to_fc (both fractal
    // coordinates). Returns undefined if the required information is not
    // yet loaded. Returns 0 if the two coordinates identify the same
    // absolute position.

    // Check for matching coordinates:
    let from_ac = fc__ac(from_fc);
    let to_ac = fc__ac(to_fc);
    if (same(from_ac, to_ac)) {
        return 0;
    }

    //  TODO: A* HERE!
    return undefined;
}


// ------------------
// Pattern Management
// ------------------

function assemble_adjacency(north, east, south, west) {
    // Takes booleans indicating north/east/south/west adjacency and
    // assembles an adjacency value (uses bit flags).
    return (
        (NORTH * north)
      | (EAST * east)
      | (SOUTH * south)
      | (WEST * west)
    );
}

function is_connected(adjacency, direction) {
    // Given an adjacency value and a direction, returns true if that
    // adjacency value includes a connection in the given direction, and
    // false otherwise.
    return !!(adjacency & direction);
}

function assemble_pattern(nw_adj, ne_adj, se_adj, sw_adj) {
    // Takes four adjacency values representing adjacency statuses of the
    // northwest, northeast, southeast, and southwest tiles of a layer,
    // and returns a single pattern value for that layer.
    return (
        nw_adj
      | (ne_adj << 4)
      | (se_adj << 8)
      | (sw_adj << 12)
    );
}

function adjacency_for(pattern, pidx) {
    // Given a pattern and a pattern index, retrieves the indicated
    // adjacency value from the pattern.
    let shift = 4 * pidx;
    return (pattern & (0xf << shift)) >>> shift;
}

function connected_in_pattern(pattern, pidx_from, pidx_to) {
    // Returns true if there is a connection from the sub-layer at the
    // given from pattern index to the sub-layer at the given to pattern
    // index within the given pattern. Only checks for a direct
    // connection, so diagonals will always return false.
    let adj = adjacency_for(pattern, pidx_from);
    if (pidx_from == 0) {
        if (pidx_to == 1) {
            return is_connected(adj, EAST);
        } else if (pidx_to == 3) {
            return is_connected(adj, SOUTH);
        } else {
            return false;
        }
    } else if (pidx_from == 1) {
        if (pidx_to == 0) {
            return is_connected(adj, WEST);
        } else if (pidx_to == 2) {
            return is_connected(adj, SOUTH);
        } else {
            return false;
        }
    } else if (pidx_from == 2) {
        if (pidx_to == 1) {
            return is_connected(adj, NORTH);
        } else if (pidx_to == 3) {
            return is_connected(adj, WEST);
        } else {
            return false;
        }
    } else if (pidx_from == 3) {
        if (pidx_to == 0) {
            return is_connected(adj, NORTH);
        } else if (pidx_to == 2) {
            return is_connected(adj, EAST);
        } else {
            return false;
        }
    } else {
        console.error("Invalid pidx", pidx_from);
    }
}

function direction_in_pattern(pattern, from_pidx, to_pidx) {
    // Returns the direction of travel necessary to take the next step
    // towards the sub-layer identified by the given "to" pattern index
    // from the sub-layer identified by the given "from" pattern index.
    // This function *assumes* that there is always full connectivity!
    // Returns 0 if the two pidx values are the same.
    if (from_pidx == to_pidx) {
        return 0;
    }

    if (connected_in_pattern(pattern, from_pidx, to_pidx)) {
        // there's a direct connection
        if (from_pidx == 0) {
            if (to_pidx == 1) { return EAST; }
            else if (to_pidx == 3) { return SOUTH; }
            else { console.error("Impossible to_pidx", to_pidx); }
        } else if (from_pidx == 1) {
            if (to_pidx == 0) { return WEST; }
            else if (to_pidx == 2) { return SOUTH; }
            else { console.error("Impossible to_pidx", to_pidx); }
        } else if (from_pidx == 2) {
            if (to_pidx == 1) { return NORTH; }
            else if (to_pidx == 3) { return WEST; }
            else { console.error("Impossible to_pidx", to_pidx); }
        } else if (from_pidx == 3) {
            if (to_pidx == 0) { return NORTH; }
            else if (to_pidx == 2) { return EAST; }
            else { console.error("Impossible to_pidx", to_pidx); }
        } else {
            console.error("Bad from_pidx", from_pidx);
        }
    } else {
        // No direct connection
        if (from_pidx == 0) {
            if (to_pidx == 1) { return SOUTH; }
            else if (to_pidx == 3) { return EAST; }
            else if (to_pidx == 2) {
                if (
                    connected_in_pattern(pattern, 0, 1)
                 && connected_in_pattern(pattern, 1, 2)
                ) { return EAST; }
                else { return SOUTH; }
            } else { console.error("Invalid to_pidx", to_pidx); }
        } else if (from_pidx == 1) {
            if (to_pidx == 0) { return SOUTH; }
            else if (to_pidx == 2) { return WEST; }
            else if (to_pidx == 3) {
                if (
                    connected_in_pattern(pattern, 1, 2)
                 && connected_in_pattern(pattern, 2, 3)
                ) { return SOUTH; }
                else { return WEST; }
            } else { console.error("Invalid to_pidx", to_pidx); }
        } else if (from_pidx == 2) {
            if (to_pidx == 1) { return WEST; }
            else if (to_pidx == 3) { return NORTH; }
            else if (to_pidx == 0) {
                if (
                    connected_in_pattern(pattern, 2, 3)
                 && connected_in_pattern(pattern, 3, 0)
                ) { return WEST; }
                else { return NORTH; }
            } else { console.error("Invalid to_pidx", to_pidx); }
        } else if (from_pidx == 3) {
            if (to_pidx == 0) { return EAST; }
            else if (to_pidx == 2) { return NORTH; }
            else if (to_pidx == 1) {
                if (
                    connected_in_pattern(pattern, 3, 0)
                 && connected_in_pattern(pattern, 0, 1)
                ) { return NORTH; }
                else { return EAST; }
            } else { console.error("Invalid to_pidx", to_pidx); }
        } else {
            console.error("Bad from_pidx", from_pidx);
        }
    }

    console.error("Fell out of cases", from_pidx, to_pidx);
}

function inner_edge(fc, edge) {
    // Returns the fractal coordinates of the layer whose edge in the
    // given direction contains the given edge and in whose parent that
    // edge is an interior edge. This will either be the specified layer
    // itself, or some ancestor, with the maximum height above the
    // current layer being proportional to the logarithm of the distance
    // from the origin (I think).

    let trace = fc_trace(fc);
    let last = trace[trace.length - 1];

    if (last == 0) {
        if (edge == EAST || edge == SOUTH) {
            return fc;
        } else {
            return inner_edge(parent_of(fc), edge);
        }
    } else if (last == 1) {
        if (edge == WEST || edge == SOUTH) {
            return fc;
        } else {
            return inner_edge(parent_of(fc), edge);
        }
    } else if (last == 2) {
        if (edge == WEST || edge == NORTH) {
            return fc;
        } else {
            return inner_edge(parent_of(fc), edge);
        }
    } else if (last == 3) {
        if (edge == EAST || edge == NORTH) {
            return fc;
        } else {
            return inner_edge(parent_of(fc), edge);
        }
    }
}

function canonical_connection_index(fc, edge) {
    // Returns the index along the indicated edge at which a connecting
    // line crosses the edge, for the canonical connection. Does not
    // determine whether there is or is not a connection crossing that
    // edge.
    let height = fc_height(fc);
    let trace = fc_trace(fc);
    let edge_height = fc_layers_below(fc);

    let edge_length = Math.pow(PATTERN_SIZE, edge_height);

    let es = edge_seed(fc, edge);

    if (edge_length == 1) {
        return 0;
    } else {
        return randint(edge_length, es);
    }
}

function all_connection_indices(fc, edge) {
    // Returns an array of connection indices, including at least the
    // canonical connection index for the given edge, and possibly other
    // indices, depending on the edge. Longer edges have a higher chance
    // of including extra connections.

    let result = [ canonical_connection_index(fc, edge) ];

    let height = fc_height(fc);
    let trace = fc_trace(fc);
    let edge_height = fc_layers_below(fc);

    let edge_length = Math.pow(PATTERN_SIZE, edge_height);

    let es = edge_seed(fc, edge);

    let max_connections = Math.max(1, Math.log2(edge_length) - 1);

    // TODO: HERE

    return result;
}


function full_adjacency_at(fc) {
    // Retrieves the full adjacency information for the given fractal
    // location, or returns undefined if there are layers that need to be
    // generated first (any required layers are enqueued for generation
    // as part of the process).
    let here_ac = fc__ac(fc);
    let here_height = fc_layers_below(fc);
    let here_edge_length = Math.pow(PATTERN_SIZE, here_height);

    let result = 0;

    for (let edge of [ NORTH, EAST, SOUTH, WEST ]) {
        let inner_fc = inner_edge(fc, edge);
        let inner_parent = parent_of(inner_fc);
        let inner_trace = fc_trace(inner_fc);
        let inner_last = inner_trace[inner_trace.length - 1];


        let parent_layer = lookup_layer(inner_parent);
        if (parent_layer == undefined) {
            return undefined;
        }

        let adj = adjacency_for(parent_layer.pattern, inner_last);

        if (is_connected(adj, edge)) {
            let cidx = canonical_connection_index(inner_fc, edge);
            let edge_ac = fc__edge_ac(inner_fc, edge);
            let ovec = ori__vec(edge_ori(edge));
            let precise_ac = [
                edge_ac[0] + cidx * ovec[0],
                edge_ac[1] + cidx * ovec[1]
            ];

            // Swap over to the correct side of NORTH and WEST edges
            if (edge == NORTH) {
                precise_ac[1] -= 1;
            } else if (edge == WEST) {
                precise_ac[0] += 1;
            }

            if (edge == NORTH || edge == SOUTH) {
                let min_x = here_ac[0];
                let max_x = here_ac[0] + here_edge_length - 1;
                if (precise_ac[0] >= min_x && precise_ac[0] <= max_x) {
                    result |= edge;
                }
            } else {
                let min_y = here_ac[1];
                let max_y = here_ac[1] + here_edge_length - 1;
                if (precise_ac[1] >= min_y && precise_ac[1] <= max_y) {
                    result |= edge;
                }
            }
        }
    }

    return result;
}

function connected_neighbors_of(seed, ac) {
    // Returns an array containing all of the absolute coordinate
    // positions which are connected and adjacent to the given absolute
    // coordinates. Returns undefined if necessary layer information has
    // not yet been generated (and the generation of that material will
    // be requested).
    let adj = full_adjacency_at(ac__fc(seed, ac));
    if (adj == undefined) {
        return undefined;
    }

    let result = [];
    if (is_connected(adj, NORTH)) {
        result.push([ac[0], ac[1] + 1]);
    }
    if (is_connected(adj, EAST)) {
        result.push([ac[0] + 1, ac[1]]);
    }
    if (is_connected(adj, SOUTH)) {
        result.push([ac[0], ac[1] - 1]);
    }
    if (is_connected(adj, WEST)) {
        result.push([ac[0] - 1, ac[1]]);
    }

    return result;
}

// The pattern for a fully-connected loop
var FULL_PATTERN = assemble_pattern(
    EAST | SOUTH,
    WEST | SOUTH,
    WEST | NORTH,
    EAST | NORTH
);
// The four possible still-connected non-loop patterns
var CUT_PATTERNS = [
    assemble_pattern(
        SOUTH,
        SOUTH,
        WEST | NORTH,
        EAST | NORTH
    ),
    assemble_pattern(
        EAST | SOUTH,
        WEST,
        WEST,
        EAST | NORTH
    ),
    assemble_pattern(
        EAST | SOUTH,
        WEST | SOUTH,
        NORTH,
        NORTH
    ),
    assemble_pattern(
        EAST,
        WEST | SOUTH,
        WEST | NORTH,
        EAST
    ),
];


// ---------
// Maze Code
// ---------

function pick_pattern(l_seed, height) {
    // Picks a pattern for a layer based on that layer's local seed.

    // Probability of cutting this loop
    let cut_prob = Math.max(0, 1.23 - 0.23 * height);
    // Pick a random pattern: either a loop, or a U-shape
    let is_cut = flip_biased(cut_prob, l_seed);

    // Pick a random cut pattern, or use the full pattern
    if (is_cut) {
        return choose_randomly(CUT_PATTERNS, l_seed);
    } else {
        return FULL_PATTERN;
    }
}

function gen_central_layer(seed, height) {
    // Using the given seed, generates and returns the central layer at
    // the given height.

    // Compute the seed
    let fc = natural_lower_coords(seed, height + 1);
    let l_seed = local_seed(fc);

    // Assemble empty result:
    let result = {
        "coords": fc,
        "seed": l_seed,
        "pattern": pick_pattern(l_seed, height),
        "children": []
    };

    // No children at height 0
    if (height == 0) {
        result.children = null;
    }

    return result;
}

function gen_layer(parent_fc, parent_layer, index) {
    // Generates a non-central layer given the fractal coordinates of its
    // parent, its parent layer, and its pattern index within that
    // parent.

    // Re-assemble full trace:
    let seed = fc_seed(parent_fc);
    let height = fc_height(parent_fc);
    let ptrace = fc_trace(parent_fc);
    let trace = ptrace.slice();
    trace.push(index);
    let fc = build_fc(seed, height, trace);
    let nr_fc = normalize_fc(fc);

    // Compute the seed
    let l_seed = local_seed(nr_fc);

    // Create result
    let result = {
        "coords": nr_fc,
        "seed": l_seed,
        "pattern": pick_pattern(l_seed, fc_layers_below(fc)),
        "children": []
    };

    if (height == 0) {
        result.children = null;
    }

    return result;
}


// -------
// Testing
// -------


var TESTS = [
    [ "blend_color:0", blend_color("#000000ff", "#ffffffff", 0.5), "#7f7f7fff"],
    [ "same:0", same([4, 4], [4, 4]), true],
    [ "same:1", same([4, 4], [5, 5]), false],
    [
        "edge_seed:0",
        edge_seed([ 17, 1, [0, 1] ], EAST),
        edge_seed([17, 1, [1, 0]], WEST)
    ],
    [
        "edge_seed:1",
        edge_seed([ 38928, 2, [2] ], NORTH),
        edge_seed([ 38928, 2, [1] ], SOUTH)
    ],
    [ "pidx__pc__pidx:0", pc__pidx(pidx__pc(0)), 0 ],
    [ "pidx__pc__pidx:1", pc__pidx(pidx__pc(1)), 1 ],
    [ "pidx__pc__pidx:2", pc__pidx(pidx__pc(2)), 2 ],
    [ "pidx__pc__pidx:3", pc__pidx(pidx__pc(3)), 3 ],
    [ "normalize_fc:0", normalize_fc([17, 2, [3, 1]]), [17, 1, [1]] ],
    [ "extend_fc:0", extend_fc([17, 1, [1, 3]]), [17, 2, [3, 1, 3]] ],
    [
        "common_parent:0",
        common_parent([17, 1, [1, 3]], [17, 1, [1, 0]]),
        [ [17, 1, [1]], 3, 0 ]
    ],
    [
        "common_parent:1",
        common_parent([17, 0, [1]], [17, 2, [0, 2, 0]]),
        [ [17, 3, [1]], 3, 0 ]
    ],
    [ "inner_edge:0", inner_edge([17, 1, [3, 3]], NORTH), [17, 1, [3, 3]] ],
    [ "inner_edge:1", inner_edge([17, 1, [3, 3]], SOUTH), [17, 3, [1]] ],
    [ "is_inside:0", is_inside([17, 3, [0, 2]], [17, 3, [0, 2, 1, 2]]), true ],
    [ "is_inside:1", is_inside([17, 3, [1, 3, 0]], [17, 1, [0, 0]]), true ],
    [ "is_inside:2", is_inside([17, 3, [1, 3, 0]], [17, 1, [0, 2]]), true ],
    [ "is_inside:3", is_inside([17, 3, [1, 3, 0]], [17, 2, [2]]), false ],
    [ "is_inside:4", is_inside([17, 1, [0]], [17, 1, [0]]), false ],
    [ "is_inside:5", is_inside([17, 1, [0]], [17, 2, [3]]), false ],
    [ "is_inside:6", is_inside([17, 1, [1]], [17, 3, [1, 3, 1, 3]]), true ],
    [
        "assemble_adjacency>-<is_connected:0",
        is_connected(assemble_adjacency(true, false, false, false), NORTH),
        true
    ],
    [
        "assemble_adjacency>-<is_connected:0",
        is_connected(assemble_adjacency(true, false, false, true), EAST),
        false
    ],
    [
        "assemble_adjacency>-<is_connected:0",
        is_connected(assemble_adjacency(true, false, false, true), SOUTH),
        false
    ],
    [
        "assemble_adjacency>-<is_connected:0",
        is_connected(assemble_adjacency(true, false, false, true), WEST),
        true
    ],
    [
        "fc__ac:0",
        fc__ac([17, 1, [2, 0]]),
        [0, -1]
    ],
    [
        "fc__ac:1",
        fc__ac([17, 1, [2, 1]]),
        [1, -1]
    ],
    [
        "fc__ac:2",
        fc__ac([17, 1, [2, 2]]),
        [1, -2]
    ],
    [
        "fc__ac:3",
        fc__ac([17, 1, [2, 3]]),
        [0, -2]
    ],
    [
        "ac__fc:0",
        ac__fc(17, [0, -1]),
        [17, 1, [2, 0]]
    ],
    [
        "ac__fc:1",
        ac__fc(17, [1, -1]),
        [17, 1, [2, 1]]
    ],
    [
        "ac__fc:2",
        ac__fc(17, [1, -2]),
        [17, 1, [2, 2]]
    ],
    [
        "ac__fc:3",
        ac__fc(17, [0, -2]),
        [17, 1, [2, 3]]
    ],
    [ "fc__ac__fc:0", fc__ac(ac__fc(17, [0, -1])), [0, -1] ],
    [ "fc__ac__fc:1", fc__ac(ac__fc(17, [1, -1])), [1, -1] ],
    [ "fc__ac__fc:2", fc__ac(ac__fc(17, [1, -2])), [1, -2] ],
    [ "fc__ac__fc:3", fc__ac(ac__fc(17, [0, -2])), [0, -2] ],
    [ "fractal_origin:0", fractal_origin(0), [0, 0] ],
    [ "fractal_origin:1", fractal_origin(1), [-2, -2] ],
    [ "fractal_origin:2", fractal_origin(2), [-2, -2] ],
    [ "fractal_origin:3", fractal_origin(3), [-10, -10] ],
    [ "fractal_origin:4", fractal_origin(4), [-10, -10] ],
    [ "fractal_origin:5", fractal_origin(5), [-42, -42] ],
    [ "fractal_origin:6", fractal_origin(6), [-42, -42] ],
    [ "fc__ac:8", fc__ac([17, 3, [2, 0]]), [-2, -6] ],
    [ "fc__ac__fc:4", fc__ac(ac__fc(17, [-2, -6])), [-2, -6] ],
    [ "fc__edge_ac:0", fc__edge_ac([17, 1, [3]], SOUTH), [-2, -2] ],
    [ "fc__edge_ac:1", fc__edge_ac([17, 3, [2, 0]], NORTH), [-2, -2] ],
    [
        "connected_in_pattern:0",
        connected_in_pattern(CUT_PATTERNS[0], 0, 1),
        false
    ],
    [
        "connected_in_pattern:1",
        connected_in_pattern(CUT_PATTERNS[0], 0, 2),
        false
    ],
    [
        "connected_in_pattern:2",
        connected_in_pattern(CUT_PATTERNS[0], 0, 3),
        true
    ],
    [
        "connected_in_pattern:3",
        connected_in_pattern(CUT_PATTERNS[0], 1, 0),
        false
    ],
    [
        "connected_in_pattern:4",
        connected_in_pattern(CUT_PATTERNS[0], 1, 2),
        true
    ],
    [
        "connected_in_pattern:5",
        connected_in_pattern(CUT_PATTERNS[0], 1, 3),
        false
    ],
    [
        "connected_in_pattern:6",
        connected_in_pattern(CUT_PATTERNS[0], 2, 0),
        false
    ],
    [
        "connected_in_pattern:7",
        connected_in_pattern(CUT_PATTERNS[0], 2, 1),
        true
    ],
    [
        "connected_in_pattern:8",
        connected_in_pattern(CUT_PATTERNS[0], 2, 3),
        true
    ],
    [
        "connected_in_pattern:9",
        connected_in_pattern(CUT_PATTERNS[0], 3, 0),
        true
    ],
    [
        "connected_in_pattern:10",
        connected_in_pattern(CUT_PATTERNS[0], 3, 1),
        false
    ],
    [
        "connected_in_pattern:11",
        connected_in_pattern(CUT_PATTERNS[0], 3, 2),
        true
    ],
    [
        "direction_in_pattern:0",
        direction_in_pattern(CUT_PATTERNS[0], 0, 1),
        SOUTH
    ],
    [
        "direction_in_pattern:1",
        direction_in_pattern(CUT_PATTERNS[0], 0, 2),
        SOUTH
    ],
    [
        "direction_in_pattern:2",
        direction_in_pattern(CUT_PATTERNS[0], 0, 3),
        SOUTH
    ],
    [
        "direction_in_pattern:3",
        direction_in_pattern(CUT_PATTERNS[0], 1, 0),
        SOUTH
    ],
    [
        "direction_in_pattern:4",
        direction_in_pattern(CUT_PATTERNS[0], 1, 3),
        SOUTH
    ],
    [
        "direction_in_pattern:5",
        direction_in_pattern(CUT_PATTERNS[0], 1, 2),
        SOUTH
    ],
    [
        "direction_in_pattern:6",
        direction_in_pattern(CUT_PATTERNS[0], 3, 0),
        NORTH
    ],
    [
        "direction_in_pattern:7",
        direction_in_pattern(CUT_PATTERNS[0], 3, 1),
        EAST
    ],
    [
        "direction_in_pattern:8",
        direction_in_pattern(CUT_PATTERNS[0], 3, 2),
        EAST
    ],
    [
        "direction_in_pattern:9",
        direction_in_pattern(CUT_PATTERNS[0], 2, 0),
        WEST
    ],
    [
        "direction_in_pattern:10",
        direction_in_pattern(CUT_PATTERNS[0], 2, 3),
        WEST
    ],
    [
        "direction_in_pattern:11",
        direction_in_pattern(CUT_PATTERNS[0], 2, 1),
        NORTH
    ],
    [
        "direction_in_pattern:12",
        direction_in_pattern(CUT_PATTERNS[1], 1, 2),
        WEST
    ],
    [
        "direction_in_pattern:13",
        direction_in_pattern(CUT_PATTERNS[1], 2, 1),
        WEST
    ],
    [
        "direction_in_pattern:14",
        direction_in_pattern(CUT_PATTERNS[1], 0, 3),
        SOUTH
    ],
    [
        "direction_in_pattern:15",
        direction_in_pattern(CUT_PATTERNS[1], 3, 0),
        NORTH
    ],
    [
        "direction_in_pattern:16",
        direction_in_pattern(CUT_PATTERNS[1], 0, 1),
        EAST
    ],
    [
        "pattern_assembly:0",
        adjacency_for(
            assemble_pattern(
                SOUTH | EAST,
                SOUTH | WEST,
                NORTH | WEST,
                NORTH | EAST
            ),
            0
        ),
        SOUTH | EAST
    ],
    [
        "pattern_assembly:1",
        adjacency_for(
            assemble_pattern(
                SOUTH | EAST,
                SOUTH | WEST,
                NORTH | WEST,
                NORTH | EAST
            ),
            1
        ),
        SOUTH | WEST
    ],
    [
        "pattern_assembly:2",
        adjacency_for(
            assemble_pattern(
                SOUTH | EAST,
                SOUTH | WEST,
                NORTH | WEST,
                NORTH | EAST
            ),
            2
        ),
        NORTH | WEST
    ],
    [
        "pattern_assembly:3",
        adjacency_for(
            assemble_pattern(
                SOUTH | EAST,
                SOUTH | WEST,
                NORTH | WEST,
                NORTH | EAST
            ),
            3
        ),
        NORTH | EAST
    ],
];


/*
 * Diagram of [X, 4, [3]] with full connectivity:
      X
    -10 -9 -8 -7 -6 -5 -4 -3 -2 -1  0  1  2  3  4  5
 
Y 5   3--3--3--3--3--3--3--3--2--2--2--2--2--2--2--2
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
  4   3--+  +--+  +--+  +--+  2--+  +--+  +--+  +--2
      |     |     |     |     |     |     |     |   
  3   3--+--+--+  +--+--+--+  2--+--+--+  +--+--+--2
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
  2   3--+  +--+  +--+  +--+  2--+  +--+  +--+  +--2
      |           |           |           |
  1   3--+--+--+--+--+--+--+  1--1--0--0--+--+--+--2
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
  0   3--+  +--+  +--+  +--+  1--+  0--0  +--+  +--2
      |     |     |     |     |     |     |     |   
 -1   3--+--+--+  +--+--+--+  1--+--+--1  +--+--+--2
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
 -2   3--+  +--+  +--+  +--+  1--1  1--1  2--2  2--2
      |                       |
 -3   3--+--+--+--+--+--+--+--+--+--+--+--+--+--+--3
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
 -4   3--+  +--+  +--+  +--+  +--+  +--+  +--+  +--3
      |     |     |     |     |     |     |     |   
 -5   3--+--+--+  +--+--+--+  +--+--+--+  +--+--+--3
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
 -6   3--+  +--+  +--+  +--+  +--+  +--+  +--+  +--3
      |           |           |           |
 -7   3--+--+--+--+--+--+--+  +--+--+--+--+--+--+--3
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
 -8   3--+  +--+  +--+  +--+  +--+  +--+  +--+  +--3
      |     |     |     |     |     |     |     |   
 -9   3--+--+--+  +--+--+--+  +--+--+--+  +--+--+--3
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
 -10  3--3  3--3  3--3  3--3  3--3  3--3  3--3  3--3
 */ 

var LATE_TESTS = [
    [
        "full_adjacency_at:0",
        [ x => is_connected(full_adjacency_at(x), SOUTH), [ [17, 0, [0]] ] ],
        [ x => is_connected(full_adjacency_at(x), NORTH), [ [17, 0, [3]] ] ],
    ],
    [
        "full_adjacency_at:1",
        [ x => is_connected(full_adjacency_at(x), EAST), [ [17, 0, [0]] ] ],
        [ x => is_connected(full_adjacency_at(x), WEST), [ [17, 0, [1]] ] ],
    ],
    [
        "full_adjacency_at:2",
        [ x => is_connected(full_adjacency_at(x), SOUTH), [ [17, 0, [1]] ] ],
        [ x => is_connected(full_adjacency_at(x), NORTH), [ [17, 0, [2]] ] ],
    ],
    [
        "full_adjacency_at:3",
        [ x => is_connected(full_adjacency_at(x), EAST), [ [17, 0, [3]] ] ],
        [ x => is_connected(full_adjacency_at(x), WEST), [ [17, 0, [2]] ] ],
    ],
    [
        "full_adjacency_at:4",
        [ x => is_connected(full_adjacency_at(x), WEST), [ [17, 0, [3]] ] ],
        [ x => is_connected(full_adjacency_at(x), EAST), [ [17, 1, [0, 1]] ] ],
    ],
    [
        "full_adjacency_at:5",
        [ x => is_connected(full_adjacency_at(x), SOUTH), [ [17, 0, [3]] ] ],
        [ x => is_connected(full_adjacency_at(x), NORTH), [ [17, 1, [2, 0]] ] ],
    ],
    [
        "full_adjacency_at:6",
        [ x => is_connected(full_adjacency_at(x), SOUTH), [ [17, 0, [2]] ] ],
        [ x => is_connected(full_adjacency_at(x), NORTH), [ [17, 1, [2, 1]] ] ],
    ],
    [
        "full_adjacency_at:7",
        [ x => is_connected(full_adjacency_at(x), EAST), [ [17, 0, [2]] ] ],
        [ x => is_connected(full_adjacency_at(x), WEST), [[17, 2, [2, 0, 3]] ]],
    ],
    [
        "full_adjacency_at:8",
        [ x => is_connected(full_adjacency_at(x), EAST), [ [17, 0, [1]] ] ],
        [ x => is_connected(full_adjacency_at(x), WEST), [[17, 2, [2, 0, 0]] ]],
    ],
    [
        "full_adjacency_at:9",
        [ x => is_connected(full_adjacency_at(x), NORTH), [ [17, 0, [1]] ] ],
        [ x => is_connected(full_adjacency_at(x), SOUTH), [[17, 2, [0, 2, 2]]]],
    ],
    [
        "full_adjacency_at:10",
        [ x => is_connected(full_adjacency_at(x), NORTH), [ [17, 0, [0]] ] ],
        [ x => is_connected(full_adjacency_at(x), SOUTH), [[17, 2, [0, 2, 3]]]],
    ],
    [
        "full_adjacency_at:11",
        [ x => is_connected(full_adjacency_at(x), WEST), [ [17, 0, [0]] ] ],
        [ x => is_connected(full_adjacency_at(x), EAST), [ [17, 1, [0, 1]] ] ],
    ],
];

function same(a, b) {
    // Object-structure based equality check
    if (Array.isArray(a)) {
        if (Array.isArray(b)) {
            if (a.length != b.length) {
                return false;
            }
            for (var i = 0; i < a.length; ++i) {
                if (!same(a[i], b[i])) {
                    return false;
                }
            }
            return true;
        } else {
            return false;
        }
    } else if (typeof a === "object") {
        if (typeof b === "object") {
            // keys & values match:
            for (var k in a) {
                if (a.hasOwnProperty(k)) {
                    if (!b.hasOwnProperty(k)) {
                        return false;
                    }
                    if (!same(a[k], b[k])) {
                        return false;
                    }
                }
            }
            // extra keys in b?
            for (var k in b) {
                if (b.hasOwnProperty(k)) {
                    if (!a.hasOwnProperty(k)) {
                        return false;
                    }
                }
            }
            return true;
        } else {
            return false;
        }
    } else {
        return a === b;
    }
}

var FAILED = false;
for (let i in TESTS) {
  let t = TESTS[i];
  let name = t[0];
  let v1 = t[1];
  let v2 = t[2];
  if (!same(v1, v2)) {
    console.error("Test '" + name + "' (#" + i + ") failed.");
    console.log("Expected:");
    console.log(v2);
    console.log("Got:");
    console.log(v1);
    FAILED = true;
  }
}

function keep_testing(tests) {
  let unresolved = [];
  for (let i in tests) {
    let t = tests[i];
    let name = t[0];
    let fv1 = t[1];
    let f1 = fv1[0];
    let a1 = fv1[1];
    let fv2 = t[2];
    let f2 = fv2[0];
    let a2 = fv2[1];

    let v1 = f1(...a1);
    let v2 = f2(...a2);

    if (v1 == undefined || v2 == undefined) {
      unresolved.push(t);
      continue;
    } else if (!same(v1, v2)) {
      console.error("Late Test '" + name + "' (#" + i + ") failed.");
      console.log("Expected:");
      console.log(v2);
      console.log("Got:");
      console.log(v1);
      FAILED = true;
    }
  }
  if (unresolved.length > 0) {
    window.setTimeout(keep_testing, TEST_DELAY, unresolved);
  } else {
    if (FAILED) {
      console.log("Late tests done failing.");
    } else {
      console.log("Late tests all passed.");
    }
  }
}

keep_testing(LATE_TESTS)

// -----
// Setup
// -----

// Run when the document is loaded unless a test failed

if (!FAILED) {
    // Grab canvas & context:
    let canvas = document.getElementById("maze");
    CTX = canvas.getContext("2d");

    // Set initial canvas size & scale:
    update_canvas_size(canvas, CTX);
    set_mode(
        document.getElementById("attract_mode").checked
        ? "attract"
        : "wait"
    );
    set_scale(CTX, 1);
    set_origin(CTX, [0, 0]);
    set_position(CTX, [0, 0]);
    set_destination(CTX, [0, 0]);
    let tempo_ctl = document.getElementById("tempo");
    set_tempo_from_slider(tempo_ctl);

    // Listen for window resizes but wait until 20 ms after the last
    // consecutive one to do anything.
    var timer_id = undefined;
    window.addEventListener("resize", function() {
        if (timer_id != undefined) {
            clearTimeout(timer_id);
            timer_id = undefined;
        }
        timer_id = setTimeout(
            function () {
                timer_id = undefined;
                update_canvas_size(canvas, CTX);
            },
            20 // milliseconds
        );
    });

    // Scrolling updates scale:
    document.onwheel = function(ev) {
        if (ev.preventDefault) { ev.preventDefault(); }
        handle_scroll(CTX, ev);
    }

    // Clicking (or tapping) sets destination:
    canvas.addEventListener("click", function (ev) { handle_tap(CTX, ev); }); 

    // Draw every frame
    window.requestAnimationFrame(draw_frame);

    // Kick off generation subsystem
    gen_step();

    // Kick off distance-measurement subsystem
    dist_step();

    // Kick off step advance
    advance_step(CTX);

    // Kick off destination scrambling
    scramble_destination(CTX);

    // Attach event handlers

    document.getElementById("tempo").addEventListener(
        "change",
        function () { set_tempo_from_slider(this); }
    );

    document.getElementById("attract_mode").addEventListener(
        "change",
        function () {
            if (this.checked) {
                set_mode("attract");
            } else {
                set_mode("wait");
            }
        }
    );
}
