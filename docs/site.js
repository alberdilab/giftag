// Page furniture shared with gifter's pkgdown site (gifter/pkgdown/extra.js),
// rewritten for Quarto's markup.

// Wordmark: "gift" in forest, "ag" muted, as gifter writes "gift" + "er".
document.addEventListener("DOMContentLoaded", () => {
  const title = document.querySelector("#quarto-header .navbar-title");

  if (!title || title.textContent.trim() !== "giftag") {
    return;
  }

  const suffix = document.createElement("span");
  suffix.className = "giftag-wordmark-suffix";
  suffix.textContent = "ag";
  title.replaceChildren("gift", suffix);
});

// Breadcrumbs: giftag > section > page, read from the navbar so they follow
// _quarto.yml. A page the navbar does not list sits directly under giftag.
document.addEventListener("DOMContentLoaded", () => {
  const main = document.querySelector("main#quarto-document-content");
  const brand = document.querySelector("#quarto-header a.navbar-brand");
  if (!main || !brand || document.body.classList.contains("giftag-home")) {
    return;
  }

  const here = (href) => {
    const url = new URL(href, window.location.href);
    return url.pathname.replace(/index\.html$/, "");
  };
  const current = here(window.location.href);
  const links = Array.from(document.querySelectorAll("#quarto-header a.nav-link, #quarto-header a.dropdown-item"))
    .filter((link) => !link.classList.contains("dropdown-toggle") && here(link.href) === current);
  const match = links.find((link) => !new URL(link.href).hash) || links[0];

  const trail = [{ label: "giftag", href: brand.href }];
  if (match) {
    const menu = match.closest(".dropdown");
    if (menu) {
      trail.push({ label: menu.querySelector(".dropdown-toggle").textContent.trim() });
    }
    trail.push({ label: match.textContent.trim() });
  } else {
    const heading = main.querySelector("h1");
    trail.push({ label: heading ? heading.textContent.trim() : document.title });
  }

  const list = document.createElement("ol");
  trail.forEach((step, index) => {
    const item = document.createElement("li");
    const last = index === trail.length - 1;
    if (step.href && !last) {
      const link = document.createElement("a");
      link.href = step.href;
      link.textContent = step.label;
      item.append(link);
    } else {
      item.textContent = step.label;
      if (last) item.setAttribute("aria-current", "page");
    }
    list.append(item);
  });

  const nav = document.createElement("nav");
  nav.className = "giftag-breadcrumbs";
  nav.setAttribute("aria-label", "Breadcrumb");
  nav.append(list);
  main.prepend(nav);
});
