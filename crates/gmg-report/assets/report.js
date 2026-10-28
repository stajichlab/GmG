(function () {
  var root = document.documentElement;
  try { var t = localStorage.getItem("gmg-theme"); if (t) root.setAttribute("data-theme", t); } catch (e) {}
  document.getElementById("theme").addEventListener("click", function () {
    var dark = root.getAttribute("data-theme") ? root.getAttribute("data-theme") === "dark" : matchMedia("(prefers-color-scheme: dark)").matches;
    var next = dark ? "light" : "dark";
    root.setAttribute("data-theme", next);
    try { localStorage.setItem("gmg-theme", next); } catch (e) {}
  });

  // tooltips: hover and keyboard focus show the same thing; every value is also in the table view
  var tip = document.getElementById("tooltip");
  function show(el, x, y) {
    var s = el.getAttribute("data-tip"); if (!s) return;
    tip.textContent = s; tip.hidden = false;
    var w = tip.offsetWidth, h = tip.offsetHeight;
    tip.style.left = Math.max(8, Math.min(x + 14, innerWidth - w - 8)) + "px";
    tip.style.top = Math.max(8, y - h - 12 < 8 ? y + 18 : y - h - 12) + "px";
  }
  document.addEventListener("mousemove", function (e) {
    var el = e.target.closest && e.target.closest("[data-tip]");
    if (el) show(el, e.clientX, e.clientY); else tip.hidden = true;
  });
  document.addEventListener("focusin", function (e) {
    var el = e.target.closest && e.target.closest("[data-tip]");
    if (el) { var r = el.getBoundingClientRect(); show(el, r.left + r.width / 2, r.top); }
  });
  document.addEventListener("focusout", function () { tip.hidden = true; });

  // table view twin for every chart
  document.querySelectorAll("[data-table-toggle]").forEach(function (b) {
    b.addEventListener("click", function () {
      var fig = b.closest(".fig"), tv = fig.querySelector(".tableview"), cb = fig.querySelector(".chartbox");
      var on = b.getAttribute("aria-pressed") !== "true";
      b.setAttribute("aria-pressed", on); b.textContent = on ? "Chart view" : "Table view";
      tv.hidden = !on; cb.hidden = on;
    });
  });

  // results table: sort + filter
  var tbody = document.querySelector("table.results tbody"), rows = Array.prototype.slice.call(tbody.rows);
  document.querySelectorAll("th[data-sort]").forEach(function (th, ci) {
    th.addEventListener("click", function () {
      var dir = th.getAttribute("data-dir") === "asc" ? "desc" : "asc";
      document.querySelectorAll("th[data-sort]").forEach(function (o) { o.removeAttribute("data-dir"); });
      th.setAttribute("data-dir", dir);
      rows.sort(function (a, b) {
        var x = a.cells[ci].getAttribute("data-v"), y = b.cells[ci].getAttribute("data-v");
        var nx = parseFloat(x), ny = parseFloat(y), c = (!isNaN(nx) && !isNaN(ny)) ? nx - ny : String(x).localeCompare(String(y));
        return dir === "asc" ? c : -c;
      });
      rows.forEach(function (r) { tbody.appendChild(r); });
    });
  });
  var fs = document.getElementById("f-status"), fq = document.getElementById("f-q"), fn = document.getElementById("f-n");
  function filt() {
    var s = fs.value, q = fq.value.trim().toLowerCase(), n = 0;
    rows.forEach(function (r) {
      var ok = (!s || r.getAttribute("data-status") === s) && (!q || r.getAttribute("data-q").indexOf(q) >= 0);
      r.hidden = !ok; if (ok) n++;
    });
    fn.textContent = n + " of " + rows.length + " inputs";
  }
  fs.addEventListener("change", filt); fq.addEventListener("input", filt); filt();

  // print: open every card so the PDF is complete
  addEventListener("beforeprint", function () { document.querySelectorAll("details.card").forEach(function (d) { d.open = true; }); });
})();
