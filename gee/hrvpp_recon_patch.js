/**
 * HR-VPP season-1 reconstruction (AMPL, LENGTH, RSLOPE, TPROD, SOS/EOS, MAXV, MINV)
 * for the area east of Sentinel-2 tile T34UGV where Copernicus HR-VPP has no data.
 *
 * Earth Engine Code Editor port of src/hrvpp_recon/gee_recon.py from
 * https://github.com/yevheniikhan/hrvpp-tile-recovery (keep both in sync).
 * Method: docs/METHOD.md; validation: docs/VALIDATION.md in that repository.
 *
 * HOW TO RUN (person with Earth Engine quota):
 *  1. Choose YOUR Cloud project in the Code Editor (top right); the exports use its quota.
 *  2. Check the settings block below (years, area, Drive folder) and press Run.
 *  3. Open the Tasks tab and press RUN for each RECON_T34UGVE_<year> task (~70 x 170 km each).
 *  4. Share the Drive folder with the person who asked for the run. Do not rename the files.
 *
 * MODE = 'test' prints values at a few pixels for a quick comparison with the Python version.
 */

// ------------------------------------------------------------------ settings
var MODE = 'export';                        // 'export' or 'test'
var YEARS = [2018, 2019, 2020, 2021, 2022, 2023, 2024];
var DRIVE_FOLDER = 'HRVPP_T34UGVE_export';
// Area in EPSG:32634 on the 10 m grid of tile T34UGV (multiples of 10 m from x=699960, y=5500020).
// Default: 10 km overlap with T34UGV, then east to ~26.1 E; 48.6-50.2 N.
var X0 = 799760, X1 = 869760, Y0 = 5390220, Y1 = 5560020;
var TILE_ID = 'T34UGVE';

// ------------------------------------------------------------------ parameters (ReconParams)
var P = {
  cs_threshold: 0.60, dvi_soil: 0.09, mdvi_cap: 0.80, leaf_g: 0.5, bin_days: 10,
  max_gap_days: 60, smooth_passes: 2, sos_frac: 0.25, eos_frac: 0.15, min_ampl: 0.05,
  min_length: 30, min_obs: 8, winter_months: [10, 11, 12, 1, 2, 3], winter_base_pct: 5
};
var FACT = {AMPL: 1e4, MAXV: 1e4, MINV: 1e4, RSLOPE: 1e4, LENGTH: 1, TPROD: 10};
var S2_SR = 'COPERNICUS/S2_SR_HARMONIZED', S2_L1C = 'COPERNICUS/S2_HARMONIZED';
var CS_PLUS = 'GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED';
var DAY = 86400000;

// ------------------------------------------------------------------ helpers (client side)
function isLeap(y) { return (y % 4 === 0 && y % 100 !== 0) || y % 400 === 0; }
function ymd(ms) { var d = new Date(ms); return d.toISOString().slice(0, 10); }
function windowBins(year) {
  var ws = Date.UTC(year - 1, 9, 1), we = Date.UTC(year + 1, 2, 31), jan1 = Date.UTC(year, 0, 1);
  var nbins = Math.floor((Math.round((we - ws) / DAY) + 1) / P.bin_days);
  var centres = [];
  for (var i = 0; i < nbins; i++) centres.push(Math.round((ws - jan1) / DAY) + P.bin_days / 2 + P.bin_days * i);
  return {ws: ws, nbins: nbins, centres: centres, jan1: jan1};
}

// ------------------------------------------------------------------ inputs
function prep(img, src) {
  var clear = img.select('cs_cdf').gte(P.cs_threshold);
  if (src === 1) {
    var scl = img.select('SCL');
    clear = clear.and(scl.neq(0)).and(scl.neq(1)).and(scl.neq(3))
                 .and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10));
  }
  var red = img.select('B4').divide(10000), nir = img.select('B8').divide(10000);
  var dvi = nir.subtract(red).rename('DVI');
  var out = dvi.updateMask(clear).addBands(ee.Image.constant(src).toByte().rename('SRC'))
               .set('sza', img.get('MEAN_SOLAR_ZENITH_ANGLE'));
  return ee.Image(out.copyProperties(img, ['system:time_start', 'system:index']));
}

function dviCollection(region, start, end) {
  var cs = ee.ImageCollection(CS_PLUS);
  var sr = ee.ImageCollection(S2_SR).filterBounds(region).filterDate(start, end)
             .linkCollection(cs, ['cs_cdf']);
  var l1c = ee.ImageCollection(S2_L1C).filterBounds(region).filterDate(start, end)
             .filter(ee.Filter.inList('system:index', sr.aggregate_array('system:index')).not())
             .linkCollection(cs, ['cs_cdf']);
  var a = sr.map(function (i) { return prep(i, 1); });
  var b = l1c.map(function (i) { return prep(i, 2); });
  return a.merge(b).filter(ee.Filter.notNull(['sza']));
}

function mdviImage(dvi, region) {
  var q95 = dvi.select('DVI').reduce(ee.Reducer.percentile([95])).rename('q95');
  var dmax = dvi.select('DVI').max().rename('dmax');
  var delta = ee.Number(dmax.subtract(q95).reduceRegion({
    reducer: ee.Reducer.mean(), geometry: region, scale: 200, maxPixels: 1e9,
    bestEffort: true, tileScale: 4}).values().get(0));
  return q95.add(delta).min(P.mdvi_cap).max(P.dvi_soil + 0.1).rename('MDVI');
}

function addPpi(dvi, mdvi) {
  var dvis = P.dvi_soil;
  return dvi.map(function (img) {
    var cosz = ee.Number(img.get('sza')).multiply(3.141592653589793 / 180).cos();
    var dc = ee.Number(0.0336).add(ee.Number(0.0477).divide(cosz));
    var qe = dc.add(ee.Number(1).subtract(dc).multiply(P.leaf_g).divide(cosz));
    var k = ee.Image(1).divide(ee.Image.constant(qe).multiply(4))
              .multiply(mdvi.add(1)).divide(ee.Image(1).subtract(mdvi));
    var d = img.select('DVI').min(mdvi.subtract(1e-4)).max(dvis);
    var ppi = k.multiply(mdvi.subtract(d).divide(mdvi.subtract(dvis)).log()).multiply(-1);
    var out = ppi.clamp(0, 3).rename('PPI').updateMask(img.select('DVI').mask()).addBands(img.select('SRC'));
    return ee.Image(out.copyProperties(img, ['system:time_start']));
  });
}

// ------------------------------------------------------------------ time series
function dekadalSeries(ppi, year) {
  var w = windowBins(year);
  var base = ppi.select('PPI').reduce(ee.Reducer.percentile([P.winter_base_pct])).rename('PPI');
  var empty = ee.Image.constant(0).toFloat().rename('PPI').updateMask(ee.Image.constant(0));
  var imgs = [];
  for (var i = 0; i < w.nbins; i++) {
    var c = w.centres[i];
    var s = w.ws + i * P.bin_days * DAY, e = s + P.bin_days * DAY;
    var centre = w.jan1 + c * DAY;
    var sub = ppi.filterDate(ymd(s), ymd(e)).select('PPI');
    var comp = ee.Image(ee.Algorithms.If(sub.size().gt(0), sub.median(), empty)).rename('PPI').toFloat();
    var nobs = ee.Image(ee.Algorithms.If(sub.size().gt(0), sub.count(), ee.Image.constant(0))).rename('N').unmask(0);
    var month = new Date(centre).getUTCMonth() + 1;
    if (P.winter_months.indexOf(month) >= 0) comp = comp.unmask(base);
    var tband = ee.Image.constant(c).toFloat().rename('T').updateMask(comp.mask());
    imgs.push(comp.addBands(tband).addBands(ee.Image(nobs).toFloat())
                  .set({'system:time_start': Math.floor(centre), 't': c, 'i': i}));
  }
  return ee.ImageCollection.fromImages(imgs);
}

function interpolateGaps(col) {
  var millis = P.max_gap_days * DAY;
  var diff = ee.Filter.maxDifference({difference: millis, leftField: 'system:time_start', rightField: 'system:time_start'});
  var afterF = ee.Filter.and(diff, ee.Filter.lessThanOrEquals({leftField: 'system:time_start', rightField: 'system:time_start'}));
  var beforeF = ee.Filter.and(diff, ee.Filter.greaterThanOrEquals({leftField: 'system:time_start', rightField: 'system:time_start'}));
  var j1 = ee.Join.saveAll({matchesKey: 'after', ordering: 'system:time_start', ascending: false})
             .apply({primary: col, secondary: col, condition: afterF});
  var j2 = ee.Join.saveAll({matchesKey: 'before', ordering: 'system:time_start', ascending: true})
             .apply({primary: j1, secondary: j1, condition: beforeF});
  return ee.ImageCollection(j2.map(function (img) {
    img = ee.Image(img);
    var bm = ee.ImageCollection.fromImages(ee.List(img.get('before'))).select(['PPI', 'T']).mosaic();
    var am = ee.ImageCollection.fromImages(ee.List(img.get('after'))).select(['PPI', 'T']).mosaic();
    var t = ee.Image.constant(img.get('t')).toFloat();
    var span = am.select('T').subtract(bm.select('T'));
    var ratio = t.subtract(bm.select('T')).divide(span.where(span.eq(0), 1));
    var interp = bm.select('PPI').add(am.select('PPI').subtract(bm.select('PPI')).multiply(ratio));
    var out = img.select('PPI').unmask(interp).addBands(img.select('N'));
    return ee.Image(out.copyProperties(img, ['system:time_start', 't', 'i']));
  })).sort('system:time_start');
}

// ------------------------------------------------------------------ metrics
function smooth121(a, n) {
  var left = a.arraySlice(0, 0, 1).arrayCat(a.arraySlice(0, 0, n - 1), 0);
  var right = a.arraySlice(0, 1, n).arrayCat(a.arraySlice(0, n - 1, n), 0);
  return left.add(a.multiply(2)).add(right).divide(4);
}

function smoothedArray(series, year) {
  var n = windowBins(year).nbins;
  var complete = series.select('PPI').count().eq(n);
  var V = series.select('PPI').map(function (i) { return i.unmask(0); }).toArray().arrayProject([0]);
  V = V.updateMask(complete);
  for (var k = 0; k < P.smooth_passes; k++) V = smooth121(V, n);
  return V;
}

function seasonMetrics(series, year) {
  var w = windowBins(year), n = w.nbins, centres = w.centres;
  var lenY = isLeap(year) ? 366 : 365, lenPrev = isLeap(year - 1) ? 366 : 365;
  var half = P.bin_days / 2;
  var V = smoothedArray(series, year);
  var nobs = series.select('N').sum().rename('NOBS');
  var idxList = []; for (var i = 0; i < n; i++) idxList.push(i);
  var idx = ee.Image(ee.Array(idxList));
  var T = ee.Image(ee.Array(centres));
  var inYear = T.gte(0).and(T.lt(lenY));
  function amin(a) { return a.arrayReduce(ee.Reducer.min(), [0]).arrayGet([0]); }
  function amax(a) { return a.arrayReduce(ee.Reducer.max(), [0]).arrayGet([0]); }
  function pick(cond, a, other) { return cond.multiply(a).add(cond.not().multiply(other)); }

  var Vm = pick(inYear, V, -1);
  var maxv = amax(Vm);
  var peak = amin(pick(Vm.eq(maxv), idx, 999)).toInt();
  var left = idx.lte(peak), right = idx.gte(peak);
  var lmin = amin(pick(left, V, 99)), rmin = amin(pick(right, V, 99));
  var lminI = amax(pick(left.and(V.eq(lmin)), idx, -1)).toInt();
  var rminI = amin(pick(right.and(V.eq(rmin)), idx, 999)).toInt();
  var minv = lmin.add(rmin).divide(2);
  var ampl = maxv.subtract(minv);
  var thrS = lmin.add(maxv.subtract(lmin).multiply(P.sos_frac));
  var thrE = rmin.add(maxv.subtract(rmin).multiply(P.eos_frac));
  var rngS = idx.gte(lminI).and(left), rngE = idx.lte(rminI).and(right);
  var sosI = amin(pick(rngS.and(V.gte(thrS)), idx, 999)).toInt();
  var eosI = amax(pick(rngE.and(V.gte(thrE)), idx, -1)).toInt();
  var s0 = sosI.subtract(1).max(lminI).clamp(0, n - 1).toInt();
  var s1 = sosI.clamp(0, n - 1).toInt();
  var e0 = eosI.clamp(0, n - 1).toInt();
  var e1 = eosI.add(1).min(rminI).clamp(0, n - 1).toInt();
  var vs0 = V.arrayGet(s0), vs1 = V.arrayGet(s1), ve0 = V.arrayGet(e0), ve1 = V.arrayGet(e1);
  var ds = vs1.subtract(vs0);
  var fs = thrS.subtract(vs0).divide(ds.where(ds.eq(0), 1)).where(ds.eq(0), 1).clamp(0, 1);
  var de = ve0.subtract(ve1);
  var fe = ve0.subtract(thrE).divide(de.where(de.eq(0), 1)).where(de.eq(0), 0).clamp(0, 1);
  var t0 = centres[0];
  var sosT = s0.add(fs.multiply(s1.subtract(s0))).multiply(P.bin_days).add(t0);
  var eosT = e0.add(fe.multiply(e1.subtract(e0))).multiply(P.bin_days).add(t0);
  var length = eosT.subtract(sosT);
  var lo = T.subtract(half).max(sosT), hi = T.add(half).min(eosT);
  var wgt = hi.subtract(lo).max(0);
  var tprod = V.multiply(wgt).arrayReduce(ee.Reducer.sum(), [0]).arrayGet([0]);
  var rslope = ve0.subtract(ve1).abs().divide(P.bin_days);
  function code(tt) {
    var tf = tt.floor();
    var yy = ee.Image(year % 100).where(tf.lt(0), (year - 1) % 100).where(tf.gte(lenY), (year + 1) % 100);
    var doy = tf.add(1).where(tf.lt(0), tf.add(lenPrev + 1)).where(tf.gte(lenY), tf.subtract(lenY - 1));
    return yy.multiply(1000).add(doy);
  }
  var valid = ampl.gte(P.min_ampl).and(length.gte(P.min_length)).and(sosI.lt(999))
      .and(eosI.gte(0)).and(sosI.lte(peak)).and(eosI.gte(peak)).and(nobs.gte(P.min_obs));
  return ee.Image.cat([ampl.rename('AMPL'), length.rename('LENGTH'), rslope.rename('RSLOPE'),
      tprod.rename('TPROD'), maxv.rename('MAXV'), minv.rename('MINV'),
      code(sosT).rename('SOSD'), code(eosT).rename('EOSD')]).updateMask(valid).addBands(nobs);
}

function reconstructYear(year, region) {
  var w = windowBins(year);
  var dvi = dviCollection(region, ymd(w.ws), ymd(Date.UTC(year + 1, 3, 1)));
  var ppi = addPpi(dvi, mdviImage(dvi, region));
  return seasonMetrics(interpolateGaps(dekadalSeries(ppi, year)), year);
}

function exportImage(m) {
  var bands = ['AMPL', 'LENGTH', 'RSLOPE', 'TPROD', 'SOSD', 'EOSD', 'MAXV', 'MINV'].map(function (b) {
    return FACT[b] ? m.select(b).multiply(FACT[b]).round() : m.select(b).round();
  });
  return ee.Image.cat(bands).addBands(m.select('NOBS')).unmask(-32768).toInt16();
}

// ------------------------------------------------------------------ run
var region = ee.Geometry.Rectangle([X0, Y0, X1, Y1], 'EPSG:32634', false);
Map.centerObject(region, 8);
Map.addLayer(region, {color: 'magenta'}, 'reconstruction area');

if (MODE === 'export') {
  var W = (X1 - X0) / 10, H = (Y1 - Y0) / 10;
  YEARS.forEach(function (y) {
    Export.image.toDrive({
      image: exportImage(reconstructYear(y, region)),
      description: 'RECON_' + TILE_ID + '_' + y,
      folder: DRIVE_FOLDER,
      fileNamePrefix: 'RECON_VPP_' + y + '_S2_' + TILE_ID + '-010m_s1',
      crs: 'EPSG:32634', crsTransform: [10, 0, X0, 0, -10, Y1],
      dimensions: W + 'x' + H, maxPixels: 1e10,
      fileDimensions: [Math.ceil(W / 256) * 256, Math.ceil(H / 256) * 256],
      formatOptions: {cloudOptimized: true}
    });
  });
  print('Tasks created for years', YEARS, 'size', W + ' x ' + H + ' px. Open the Tasks tab and press RUN.');
} else {
  // Test: 2021, MDVI over tile T34UGV (as in the Python run), values at 6 pixels of the overlap with T34UFV.
  var tile = ee.Geometry.Rectangle([699960, 5390220, 809760, 5500020], 'EPSG:32634', false);
  var w21 = windowBins(2021);
  var dvi = dviCollection(tile, ymd(w21.ws), ymd(Date.UTC(2022, 3, 1)));
  var m = exportImage(seasonMetrics(interpolateGaps(dekadalSeries(addPpi(dvi, mdviImage(dvi, tile)), 2021)), 2021));
  var pts = ee.FeatureCollection([[702105, 5451895], [703505, 5450505], [704905, 5449105],
                                   [702805, 5448405], [705605, 5450905], [703005, 5452005 - 1300]]
    .map(function (xy, i) { return ee.Feature(ee.Geometry.Point(xy, 'EPSG:32634'), {id: i}); }));
  print('n scenes', dvi.size());
  print(m.sampleRegions({collection: pts, scale: 10, projection: 'EPSG:32634', tileScale: 4}));
}
