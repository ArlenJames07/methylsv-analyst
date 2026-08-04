const formatFileSize = (bytes) => {
    if (!Number.isFinite(bytes) || bytes < 0) {
        return "";
    }

    if (bytes < 1024) {
        return `${bytes} B`;
    }

    const units = ["KB", "MB", "GB", "TB"];
    let size = bytes / 1024;
    let unitIndex = 0;

    while (size >= 1024 && unitIndex < units.length - 1) {
        size /= 1024;
        unitIndex += 1;
    }

    const precision = size >= 10 ? 0 : 1;
    return `${size.toFixed(precision)} ${units[unitIndex]}`;
};

document.querySelectorAll("[data-file-input]").forEach((input) => {
    const control = input.closest(".file-control");
    const label = control?.querySelector("[data-file-label]");
    const hint = control?.querySelector(".file-copy small");
    const initialLabel = label?.textContent;
    const initialHint = hint?.textContent;

    input.addEventListener("change", () => {
        const [file] = input.files;

        if (label) {
            label.textContent = file?.name || initialLabel;
        }

        if (hint) {
            hint.textContent = file
                ? `${formatFileSize(file.size)} · ready to validate`
                : initialHint;
        }
    });
});

document.querySelectorAll("[data-upload-form]").forEach((form) => {
    form.addEventListener("submit", () => {
        const button = form.querySelector("[data-submit-button]");
        const text = form.querySelector("[data-submit-text]");

        if (!button || !form.checkValidity()) {
            return;
        }

        button.setAttribute("aria-busy", "true");
        button.disabled = true;

        if (text) {
            text.textContent = "Validating files…";
        }
    });
});
